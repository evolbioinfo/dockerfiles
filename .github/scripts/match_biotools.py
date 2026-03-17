#!/usr/bin/env python3
"""match_biotools.py – Match tools in the evolbioinfo/dockerfiles repository
with their identifiers in the bio.tools database.

The script discovers all tool directories that contain at least one Dockerfile,
queries the bio.tools REST API (https://bio.tools) for each tool, and writes a
mapping of directory names to bio.tools identifiers and descriptions.

Usage
-----
    python .github/scripts/match_biotools.py [options]

Options
-------
    --repo-dir DIR      Path to the repository root (default: three levels above
                        this script, i.e. the repository root).
    --output FILE       Write results to FILE instead of stdout.
    --format {csv,tsv,json,markdown}
                        Output format (default: csv).
    --delay SECONDS     Seconds to wait between API requests (default: 0.5).
    --no-cache          Do not cache API responses on disk.
    --cache-dir DIR     Directory used for the optional on-disk cache
                        (default: <repo-dir>/.biotools_cache).
    --list-tools        Print discovered tool directories and exit without
                        querying bio.tools.

Examples
--------
    # Print CSV results to stdout
    python .github/scripts/match_biotools.py

    # Save Markdown output to BIOTOOLS.md
    python .github/scripts/match_biotools.py --format markdown --output BIOTOOLS.md

    # Save JSON output to a file
    python .github/scripts/match_biotools.py --format json --output biotools_matches.json

    # Just list the tools that would be processed
    python .github/scripts/match_biotools.py --list-tools
"""

import argparse
import csv
import json
import os
import sys
import time
from datetime import datetime, timezone
from io import StringIO
from pathlib import Path

try:
    import requests
except ImportError:
    print(
        "Error: the 'requests' package is required.\n"
        "Install it with:  pip install requests",
        file=sys.stderr,
    )
    sys.exit(1)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

BIOTOOLS_API = "https://bio.tools/api"
BIOTOOLS_BASE = "https://bio.tools"
REPO_URL = "https://github.com/evolbioinfo/dockerfiles"

# Directories that are *not* stand-alone bioinformatics tools (base images,
# language environments, generic utilities) and that are unlikely to have a
# bio.tools entry.  They are still processed when --include-non-tools is given.
NON_TOOL_DIRS = {
    "jq",
    "perl",
    "python",
    "python-dl",
    "python-evol",
    "python-ml",
    "r-base",
    "r-evol",
    "r-extended",
    "r-gisaid",
    "r-sra",
    "s3utils",
    "sphinx",
    "ubuntu",
    "wget",
}

# ---------------------------------------------------------------------------
# Tool discovery
# ---------------------------------------------------------------------------


def discover_tools(base_dir: Path) -> list[str]:
    """Return a sorted list of tool directory names that contain a Dockerfile."""
    tools = []
    for entry in sorted(base_dir.iterdir()):
        if not entry.is_dir() or entry.name.startswith("."):
            continue
        if any(entry.rglob("Dockerfile")):
            tools.append(entry.name)
    return tools


# ---------------------------------------------------------------------------
# Bio.tools API helpers
# ---------------------------------------------------------------------------


def _get(url: str, params: dict, timeout: int = 10) -> dict | None:
    """Perform a GET request and return parsed JSON, or None on failure."""
    try:
        response = requests.get(url, params=params, timeout=timeout)
        if response.status_code == 200:
            return response.json()
        if response.status_code == 404:
            return None
        response.raise_for_status()
    except requests.RequestException as exc:
        print(f"  [warning] Request failed: {exc}", file=sys.stderr)
    return None


def _direct_lookup(tool_id: str) -> dict | None:
    """Fetch a tool entry directly by its bio.tools ID."""
    return _get(f"{BIOTOOLS_API}/tool/{tool_id}/", {"format": "json"})


def _name_search(name: str) -> dict | None:
    """Search bio.tools by exact tool name."""
    return _get(f"{BIOTOOLS_API}/tool/", {"name": name, "format": "json"})


def _text_search(query: str) -> dict | None:
    """Search bio.tools using a free-text query."""
    return _get(f"{BIOTOOLS_API}/tool/", {"q": query, "format": "json"})


# ---------------------------------------------------------------------------
# Name normalization variants
# ---------------------------------------------------------------------------


def _name_variants(name: str) -> list[str]:
    """Return a list of alternative spellings to try for a tool name."""
    variants: list[str] = [name]
    # underscore ↔ hyphen
    if "_" in name:
        variants.append(name.replace("_", "-"))
    if "-" in name:
        variants.append(name.replace("-", "_"))
    # drop separators entirely
    clean = name.replace("_", "").replace("-", "")
    if clean not in variants:
        variants.append(clean)
    return variants


# ---------------------------------------------------------------------------
# Text helpers
# ---------------------------------------------------------------------------


def _clean_description(text: str) -> str:
    """Replace newline characters with a single space and strip leading/trailing whitespace."""
    return " ".join(text.splitlines()).strip()


# ---------------------------------------------------------------------------
# Matching logic
# ---------------------------------------------------------------------------


def _pick_best(results: dict, candidates: list[str]) -> tuple[str, str, str, str] | None:
    """
    From a search result payload, pick the best matching entry.

    Returns (biotoolsID, name, description, match_type) or None.
    Priority:
      1. biotoolsID exactly matches one of *candidates*
      2. Tool name (case-insensitive) matches one of *candidates*
      3. Single result in the response list
    """
    tool_list = results.get("list") or []
    if not tool_list:
        return None

    candidates_lower = [c.lower() for c in candidates]

    # Priority 1 – biotoolsID exact match
    for t in tool_list:
        if t.get("biotoolsID", "").lower() in candidates_lower:
            return t["biotoolsID"], t.get("name", ""), _clean_description(t.get("description", "")), "exact_id"

    # Priority 2 – tool name exact match
    for t in tool_list:
        if t.get("name", "").lower() in candidates_lower:
            return t["biotoolsID"], t.get("name", ""), _clean_description(t.get("description", "")), "exact_name"

    # Priority 3 – single result (likely match)
    if results.get("count") == 1:
        t = tool_list[0]
        return t["biotoolsID"], t.get("name", ""), _clean_description(t.get("description", "")), "single_result"

    return None


def find_match(
    tool_dir: str,
    delay: float,
    cache: dict,
) -> dict:
    """
    Return a result dict with keys: tool, biotools_id, biotools_name,
    description, match_type.

    Matching strategy (stops at first success):
      1. Direct lookup  GET /api/tool/{tool_dir}/
      2. Direct lookup  GET /api/tool/{variant}/   for each name variant
      3. Name search    GET /api/tool/?name={variant}  for each variant
      4. Free-text      GET /api/tool/?q={tool_dir}
    """
    result_template = {
        "tool": tool_dir,
        "biotools_id": "",
        "biotools_name": "",
        "description": "",
        "match_type": "not_found",
    }

    if tool_dir in cache:
        return cache[tool_dir]

    variants = _name_variants(tool_dir)

    # --- Strategy 1 & 2: direct ID lookup ---
    for variant in variants:
        data = _direct_lookup(variant)
        time.sleep(delay)
        if data and data.get("biotoolsID"):
            result = {
                **result_template,
                "biotools_id": data["biotoolsID"],
                "biotools_name": data.get("name", ""),
                "description": _clean_description(data.get("description", "")),
                "match_type": "direct",
            }
            cache[tool_dir] = result
            return result

    # --- Strategy 3: name search ---
    for variant in variants:
        data = _name_search(variant)
        time.sleep(delay)
        if data:
            match = _pick_best(data, variants)
            if match:
                biotools_id, biotools_name, description, match_type = match
                result = {
                    **result_template,
                    "biotools_id": biotools_id,
                    "biotools_name": biotools_name,
                    "description": description,
                    "match_type": match_type,
                }
                cache[tool_dir] = result
                return result

    # --- Strategy 4: free-text search ---
    data = _text_search(tool_dir)
    time.sleep(delay)
    if data:
        match = _pick_best(data, variants)
        if match:
            biotools_id, biotools_name, description, match_type = match
            result = {
                **result_template,
                "biotools_id": biotools_id,
                "biotools_name": biotools_name,
                "description": description,
                "match_type": "text_search",
            }
            cache[tool_dir] = result
            return result

    cache[tool_dir] = result_template
    return result_template


# ---------------------------------------------------------------------------
# Disk cache helpers
# ---------------------------------------------------------------------------


def load_cache(cache_dir: Path) -> dict:
    cache_file = cache_dir / "results.json"
    if cache_file.exists():
        try:
            return json.loads(cache_file.read_text())
        except (json.JSONDecodeError, OSError):
            pass
    return {}


def save_cache(cache_dir: Path, cache: dict) -> None:
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_file = cache_dir / "results.json"
    cache_file.write_text(json.dumps(cache, indent=2))


# ---------------------------------------------------------------------------
# Output formatters
# ---------------------------------------------------------------------------

FIELDNAMES = ["tool", "biotools_id", "biotools_name", "description", "match_type"]


def format_csv(results: list[dict], delimiter: str = ",") -> str:
    buf = StringIO()
    writer = csv.DictWriter(buf, fieldnames=FIELDNAMES, delimiter=delimiter)
    writer.writeheader()
    writer.writerows(results)
    return buf.getvalue()


def format_json(results: list[dict]) -> str:
    return json.dumps(results, indent=2)


def _md_escape(text: str) -> str:
    """Escape pipe characters so they don't break Markdown table cells."""
    return text.replace("|", "\\|")


def format_markdown(results: list[dict]) -> str:
    """
    Render results as a Markdown document suitable for committing as BIOTOOLS.md.

    The document contains:
      - A title and short description
      - The generation timestamp
      - A summary line (matched / total)
      - A Markdown table with one row per tool
    """
    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    found = sum(1 for r in results if r["biotools_id"])
    total = len(results)

    lines: list[str] = [
        "# bio.tools Mapping",
        "",
        "This table maps the bioinformatics tools available in the "
        f"[evolbioinfo/dockerfiles]({REPO_URL}) repository to their entries in "
        f"the [bio.tools](https://bio.tools) registry.",
        "",
        f"**Generated:** {now}  ",
        f"**Matched:** {found} / {total} tools",
        "",
        "| Tool | bio.tools ID | Name | Description |",
        "|------|-------------|------|-------------|",
    ]

    for r in results:
        tool = r["tool"]
        biotools_id = r["biotools_id"]
        biotools_name = _md_escape(r["biotools_name"])
        description = _md_escape(r["description"])

        tool_cell = f"[{tool}]({REPO_URL}/tree/main/{tool})"
        if biotools_id:
            id_cell = f"[{biotools_id}]({BIOTOOLS_BASE}/{biotools_id})"
        else:
            id_cell = ""

        lines.append(f"| {tool_cell} | {id_cell} | {biotools_name} | {description} |")

    lines.append("")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Match tools in the evolbioinfo/dockerfiles repository with their "
            "identifiers in the bio.tools database."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__.split("Usage")[0],
    )
    parser.add_argument(
        "--repo-dir",
        default=None,
        metavar="DIR",
        help=(
            "Path to the repository root. Defaults to three levels above this "
            "script (i.e. the repository root when the script lives in "
            ".github/scripts/)."
        ),
    )
    parser.add_argument(
        "--output",
        default=None,
        metavar="FILE",
        help="Write results to FILE instead of stdout.",
    )
    parser.add_argument(
        "--format",
        choices=["csv", "tsv", "json", "markdown"],
        default="csv",
        help="Output format (default: csv).",
    )
    parser.add_argument(
        "--delay",
        type=float,
        default=0.5,
        metavar="SECONDS",
        help="Seconds to wait between API requests (default: 0.5).",
    )
    parser.add_argument(
        "--no-cache",
        action="store_true",
        help="Do not read from or write to the on-disk cache.",
    )
    parser.add_argument(
        "--cache-dir",
        default=None,
        metavar="DIR",
        help=(
            "Directory used for the optional on-disk response cache "
            "(default: <repo-dir>/.biotools_cache)."
        ),
    )
    parser.add_argument(
        "--include-non-tools",
        action="store_true",
        help=(
            "Include base-image and utility directories (e.g. ubuntu, python, "
            "jq) in the bio.tools search."
        ),
    )
    parser.add_argument(
        "--list-tools",
        action="store_true",
        help="Print discovered tool directories and exit without querying bio.tools.",
    )
    return parser


def main(argv: list[str] | None = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)

    # Resolve repository directory
    # Default: three levels up from this script (.github/scripts/ → repo root)
    if args.repo_dir:
        repo_dir = Path(args.repo_dir)
    else:
        repo_dir = Path(__file__).parent.parent.parent
    if not repo_dir.is_dir():
        print(f"Error: repository directory not found: {repo_dir}", file=sys.stderr)
        sys.exit(1)

    # Discover tool directories
    all_tools = discover_tools(repo_dir)
    if not args.include_non_tools:
        tools = [t for t in all_tools if t not in NON_TOOL_DIRS]
    else:
        tools = all_tools

    if args.list_tools:
        for t in tools:
            print(t)
        return

    # Cache setup
    use_cache = not args.no_cache
    cache_dir = Path(args.cache_dir) if args.cache_dir else repo_dir / ".biotools_cache"
    cache: dict = load_cache(cache_dir) if use_cache else {}

    print(
        f"Processing {len(tools)} tool directories against bio.tools …",
        file=sys.stderr,
    )

    results = []
    for idx, tool in enumerate(tools, start=1):
        if tool in cache:
            print(f"  [{idx:3d}/{len(tools)}] {tool} (cached)", file=sys.stderr)
            result = cache[tool]
        else:
            print(f"  [{idx:3d}/{len(tools)}] {tool}", file=sys.stderr)
            result = find_match(tool, delay=args.delay, cache=cache)
            # Persist cache after every new lookup so partial results are not lost
            if use_cache:
                save_cache(cache_dir, cache)
        results.append(result)

    # Produce output
    if args.format == "json":
        output = format_json(results)
    elif args.format == "tsv":
        output = format_csv(results, delimiter="\t")
    elif args.format == "markdown":
        output = format_markdown(results)
    else:
        output = format_csv(results)

    if args.output:
        Path(args.output).write_text(output)
        print(f"Results written to {args.output}", file=sys.stderr)
    else:
        print(output)

    # Summary
    found = sum(1 for r in results if r["biotools_id"])
    match_type_counts: dict[str, int] = {}
    for r in results:
        mt = r["match_type"]
        match_type_counts[mt] = match_type_counts.get(mt, 0) + 1

    print(
        f"\nSummary: {found}/{len(results)} tools matched in bio.tools",
        file=sys.stderr,
    )
    for mt, count in sorted(match_type_counts.items()):
        print(f"  {mt}: {count}", file=sys.stderr)


if __name__ == "__main__":
    main()
