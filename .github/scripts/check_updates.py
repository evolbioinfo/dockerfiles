#!/usr/bin/env python3
"""
Check for newer versions of tools in the evolbioinfo/dockerfiles repository.

For each tool, this script:
1. Finds the latest version directory in the tool's folder
2. Reads the Dockerfile to detect the upstream source (GitHub, GitLab, PyPI, Conda)
3. Queries the upstream API for the latest version
4. Compares to the current version and reports outdated tools
"""

import os
import re
import sys
import json
import urllib.request
import urllib.error
from pathlib import Path
from packaging.version import Version, InvalidVersion


REPO_ROOT = Path(__file__).parent.parent.parent

# Directories to skip (not real tools, or version-based directories)
SKIP_DIRS = {".git", ".github", "README.md"}

# Directories that use non-release versioning (commits, custom schemes)
SKIP_VERSION_CHECK = {
    "artic-ncov2019",  # commit hash
    "ml_bootstrap",  # commit hash
    "ngphylogeny_multitools",  # branch name
    "pcoc",  # commit hash
    "polecat",  # commit hash
    "reseq",  # commit hash
    "snag",  # 'master'
    "strainline",  # commit hash
    "table2itol",  # commit hash
    "treedater",  # commit hash
    "treestructure",  # commit hash
    "inkscape",  # 'latest' tag
    "sdrmhunter",  # custom script; pip deps are not the tool itself
    # Base / environment images that bundle many packages – no single tool version
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
    "ubuntu",
}

# Generic utility/dependency packages that are commonly installed alongside
# the main tool but are not the tool itself.  Excluded from pip & conda
# package detection so we pick the right primary package.
COMMON_UTILITY_PACKAGES = {
    # pip
    "pip", "setuptools", "wheel", "cryptography",
    "numpy", "scipy", "pandas", "matplotlib", "biopython",
    "pysam", "docutils", "gql", "ete3", "requests",
    "cachetools", "decorator",
    # conda
    "python", "mamba", "conda",
    "prodigal", "hmmer", "pplacer",
}

# Preferred order of conda channels for version resolution.
# Earlier entries are preferred over later ones.
CONDA_CHANNEL_PRIORITY = ["bioconda", "conda-forge", "anaconda", "defaults"]


def get_tool_dirs():
    """Return a sorted list of tool directory names."""
    dirs = []
    for entry in REPO_ROOT.iterdir():
        if entry.is_dir() and entry.name not in SKIP_DIRS and not entry.name.startswith("."):
            dirs.append(entry.name)
    return sorted(dirs)


def _version_sort_key(version_str):
    """
    Return a sort key for a version string, handling non-semver suffixes.

    Strips architecture/build suffixes (e.g. -sse3, -avx2) and then
    parses via packaging.version.Version. Falls back to the original string.
    Commit hashes (pure hex strings) are sorted before real versions.
    """
    # Identify pure git commit hashes (7-40 hex characters)
    if re.fullmatch(r"[0-9a-f]{7,40}", version_str, re.IGNORECASE):
        return (0, version_str)
    # Strip leading v/V
    v = version_str.lstrip("vV")
    # Remove known build suffixes that are not part of the semantic version
    v = re.sub(r"[-_](sse\d*|avx\d*|pthreads?|mpi|gcc|clang|linux|osx|win).*$", "", v, flags=re.IGNORECASE)
    try:
        return (1, Version(v))
    except InvalidVersion:
        return (0, version_str)


def get_latest_version_dir(tool_dir):
    """
    Return the latest versioned subdirectory for a tool.
    Skips 'dev' and 'latest' symlinks/directories.
    """
    version_dirs = []
    for entry in tool_dir.iterdir():
        if entry.is_dir() and entry.name not in ("dev", "latest"):
            version_dirs.append(entry.name)
    if not version_dirs:
        return None
    version_dirs.sort(key=_version_sort_key)
    return version_dirs[-1]


def extract_github_url(dockerfile_path, tool_name=None):
    """
    Extract the first GitHub repository URL from a Dockerfile.

    Looks for:
    - Header comment: # https://github.com/owner/repo
    - git clone https://github.com/owner/repo
    - wget .../github.com/owner/repo/...

    If tool_name is provided, prefers a match where the repo name contains
    the tool name (useful when multiple GitHub URLs appear in one Dockerfile).
    """
    try:
        content = dockerfile_path.read_text(errors="replace")
    except OSError:
        return None

    # Pattern to extract owner/repo from any github.com URL
    # Matches https://github.com/owner/repo (with optional trailing path)
    pattern = re.compile(
        r"https://github\.com/([A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+)"
    )

    def clean_repo_path(raw):
        """Return 'owner/repo' stripping any trailing .git or extra path."""
        parts = raw.rstrip("/").split("/")
        owner = parts[0]
        repo = parts[1].removesuffix(".git")
        return f"{owner}/{repo}"

    # First pass: collect all GitHub repos from header comments
    header_repos = []
    for line in content.splitlines():
        if line.startswith("#"):
            match = pattern.search(line)
            if match:
                header_repos.append(clean_repo_path(match.group(1)))

    if header_repos:
        # If there's a tool name hint, prefer repos whose name matches
        if tool_name:
            tool_key = tool_name.lower().replace("-", "").replace("_", "")
            for repo in header_repos:
                repo_name = repo.split("/")[1].lower().replace("-", "").replace("_", "")
                if tool_key in repo_name or repo_name in tool_key:
                    return repo
        return header_repos[0]

    # Second pass: collect all GitHub repos from the full file body
    # Only return a match if we can confidently link it to the tool
    body_repos = []
    for line in content.splitlines():
        match = pattern.search(line)
        if match:
            body_repos.append(clean_repo_path(match.group(1)))

    if body_repos and tool_name:
        tool_key = tool_name.lower().replace("-", "").replace("_", "")
        for repo in body_repos:
            repo_name = repo.split("/")[1].lower().replace("-", "").replace("_", "")
            if tool_key in repo_name or repo_name in tool_key:
                return repo

    # No reliable match found
    return None


def extract_gitlab_url(dockerfile_path, tool_name=None):
    """
    Extract a GitLab repository URL from a Dockerfile (gitlab.com hosted).

    Returns 'namespace/project' suitable for the GitLab API, or None.
    """
    try:
        content = dockerfile_path.read_text(errors="replace")
    except OSError:
        return None

    pattern = re.compile(
        r"https://gitlab\.com/([A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+)"
    )

    def clean_path(raw):
        parts = raw.rstrip("/").split("/")
        return f"{parts[0]}/{parts[1]}"

    candidates = []
    for line in content.splitlines():
        match = pattern.search(line)
        if match:
            candidates.append(clean_path(match.group(1)))

    if not candidates:
        return None

    if tool_name:
        tool_key = tool_name.lower().replace("-", "").replace("_", "")
        for repo in candidates:
            repo_name = repo.split("/")[1].lower().replace("-", "").replace("_", "")
            if tool_key in repo_name or repo_name in tool_key:
                return repo

    return candidates[0]


def extract_pypi_package(dockerfile_path, tool_name=None):
    """
    Extract the PyPI package name from a `pip install <pkg>==<version>` line.

    Returns the package name (str) that most closely matches tool_name,
    or None if no versioned pip install is found.
    """
    try:
        content = dockerfile_path.read_text(errors="replace")
    except OSError:
        return None

    # Match: pip[3] install [opts] pkg==version [more pkgs...]
    pkg_pattern = re.compile(r"([A-Za-z0-9_.-]+)==([0-9][A-Za-z0-9._-]*)")
    pip_line_pattern = re.compile(r"pip3?\s+install\b")

    packages = []
    for line in content.splitlines():
        if pip_line_pattern.search(line):
            for m in pkg_pattern.finditer(line):
                pkg = m.group(1)
                if pkg.lower() not in COMMON_UTILITY_PACKAGES:
                    packages.append(pkg)

    if not packages:
        return None

    if tool_name:
        tool_key = tool_name.lower().replace("-", "").replace("_", "")
        # Exact or substring match first
        for pkg in packages:
            pkg_key = pkg.lower().replace("-", "").replace("_", "")
            if tool_key == pkg_key or tool_key in pkg_key or pkg_key in tool_key:
                return pkg
    # Fallback: return the first (and often only) package
    return packages[0]


def extract_conda_package(dockerfile_path, tool_name=None):
    """
    Extract the conda/bioconda package name from a conda/mamba install line.

    Returns (package_name, channel) or None.
    Channels are resolved using CONDA_CHANNEL_PRIORITY.
    """
    try:
        content = dockerfile_path.read_text(errors="replace")
    except OSError:
        return None

    # Match: mamba/conda install [-c channel] pkg=version [more...]
    install_pattern = re.compile(
        r"(?:mamba|conda)\s+install\b(.*?)(?=\\|$)", re.IGNORECASE
    )
    pkg_pattern = re.compile(r"([A-Za-z0-9_.-]+)=([0-9][A-Za-z0-9._-]*)")
    channel_pattern = re.compile(r"-c\s+(\S+)")

    candidates = []  # list of (pkg, channel)
    for line in content.splitlines():
        m = install_pattern.search(line)
        if not m:
            continue
        args = m.group(1)
        channels = channel_pattern.findall(args)
        channel = "defaults"
        for preferred in CONDA_CHANNEL_PRIORITY:
            if preferred in channels:
                channel = preferred
                break
        else:
            if channels:
                channel = channels[0]
        for pm in pkg_pattern.finditer(args):
            pkg = pm.group(1)
            if pkg.lower() not in COMMON_UTILITY_PACKAGES:
                candidates.append((pkg, channel))

    if not candidates:
        return None

    if tool_name:
        tool_key = tool_name.lower().replace("-", "").replace("_", "")
        for pkg, channel in candidates:
            pkg_key = pkg.lower().replace("-", "").replace("_", "")
            if tool_key == pkg_key or tool_key in pkg_key or pkg_key in tool_key:
                return (pkg, channel)

    return candidates[0]


def http_get_json(url, headers=None):
    """Make an HTTP GET request and return the parsed JSON body, or None."""
    req = urllib.request.Request(url)
    req.add_header("User-Agent", "evolbioinfo-dockerfiles-checker")
    if headers:
        for k, v in headers.items():
            req.add_header(k, v)
    try:
        with urllib.request.urlopen(req, timeout=15) as resp:
            return json.loads(resp.read().decode())
    except urllib.error.HTTPError as e:
        if e.code == 404:
            return None
        if e.code in (403, 429):
            print(f"Warning: rate limit or auth error for {url} (HTTP {e.code})", file=sys.stderr)
            return None
        print(f"Warning: HTTP {e.code} for {url}", file=sys.stderr)
        return None
    except (urllib.error.URLError, OSError) as e:
        print(f"Warning: Network error for {url}: {e}", file=sys.stderr)
        return None


def get_pypi_latest_version(package_name):
    """Return the latest version of a PyPI package, or None."""
    url = f"https://pypi.org/pypi/{package_name}/json"
    data = http_get_json(url)
    if data and "info" in data:
        return data["info"]["version"]
    return None


def get_conda_latest_version(package_name, channel="bioconda"):
    """
    Return the latest version of a package in a conda channel, or None.

    Queries the Anaconda.org API.
    """
    url = f"https://api.anaconda.org/package/{channel}/{package_name}"
    data = http_get_json(url)
    if data and "latest_version" in data:
        return data["latest_version"]
    return None


def get_gitlab_latest_version(namespace_project):
    """
    Return the latest tag for a gitlab.com project, or None.

    namespace_project is 'namespace/project' (e.g. 'sysimm/mafft').
    """
    encoded = namespace_project.replace("/", "%2F")
    url = f"https://gitlab.com/api/v4/projects/{encoded}/releases?per_page=1"
    data = http_get_json(url)
    if data and isinstance(data, list) and data:
        return data[0].get("tag_name")
    # Fall back to tags API
    url = f"https://gitlab.com/api/v4/projects/{encoded}/repository/tags?per_page=1"
    data = http_get_json(url)
    if data and isinstance(data, list) and data:
        return data[0].get("name")
    return None


def github_api_request(url, token=None):
    """Make a GitHub API request and return the parsed JSON."""
    req = urllib.request.Request(url)
    req.add_header("Accept", "application/vnd.github+json")
    req.add_header("X-GitHub-Api-Version", "2022-11-28")
    req.add_header("User-Agent", "evolbioinfo-dockerfiles-checker")
    if token:
        req.add_header("Authorization", f"Bearer {token}")
    try:
        with urllib.request.urlopen(req, timeout=15) as resp:
            return json.loads(resp.read().decode())
    except urllib.error.HTTPError as e:
        if e.code in (403, 429):
            # Rate limited or forbidden without a token
            print(
                f"Warning: GitHub API rate limit or auth error for {url} "
                f"(HTTP {e.code}). Set GITHUB_TOKEN to avoid rate limiting.",
                file=sys.stderr,
            )
            return None
        if e.code == 404:
            return None
        print(f"Warning: GitHub API error {e.code} for {url}", file=sys.stderr)
        return None
    except (urllib.error.URLError, OSError) as e:
        print(f"Warning: Network error for {url}: {e}", file=sys.stderr)
        return None


def get_latest_release(repo, token=None):
    """Return the tag name of the latest GitHub release, or None."""
    url = f"https://api.github.com/repos/{repo}/releases/latest"
    data = github_api_request(url, token)
    if data and "tag_name" in data:
        return data["tag_name"]
    return None


def get_latest_tag(repo, token=None):
    """Return the name of the latest git tag, or None."""
    url = f"https://api.github.com/repos/{repo}/tags?per_page=1"
    data = github_api_request(url, token)
    if data and isinstance(data, list) and data:
        return data[0]["name"]
    return None


def normalize_version(version_str):
    """
    Normalize a version string for comparison.
    Returns a Version object, or None if it cannot be parsed.
    """
    if not version_str:
        return None
    # Strip common prefixes
    v = version_str.lstrip("vV")
    # Remove build metadata or commit suffixes (e.g. v2.3.3-1-gabcdef)
    v = re.sub(r"-\d+-g[0-9a-f]+$", "", v)
    try:
        return Version(v)
    except InvalidVersion:
        return None


def is_newer(upstream_tag, current_dir_name):
    """
    Return True if the upstream_tag is newer than current_dir_name.
    Falls back to string comparison if version parsing fails.
    """
    upstream_v = normalize_version(upstream_tag)
    # The current version is encoded in the directory name
    current_v = normalize_version(current_dir_name)

    if upstream_v is not None and current_v is not None:
        return upstream_v > current_v

    # Cannot parse versions – do a plain string comparison
    upstream_clean = upstream_tag.lstrip("vV")
    current_clean = current_dir_name.lstrip("vV")
    return upstream_clean != current_clean


def check_tool(tool_name, token=None):
    """
    Check a single tool for available updates.

    Tries upstream sources in order:
    1. GitHub (releases → tags)
    2. GitLab (gitlab.com, releases → tags)
    3. PyPI (pip install <pkg>==<version>)
    4. Conda / Bioconda (conda/mamba install <pkg>=<version>)

    Returns a dict with:
      - tool: tool name
      - current: current version directory name
      - upstream: latest upstream version string
      - source: where the upstream version was found ('github', 'gitlab', 'pypi', 'conda')
      - source_ref: human-readable reference (repo path, package name, etc.)
      - outdated: True if upstream > current
      - error: error message if something went wrong (or None)
    """
    if tool_name in SKIP_VERSION_CHECK:
        return None

    tool_dir = REPO_ROOT / tool_name
    if not tool_dir.is_dir():
        return None

    current_version = get_latest_version_dir(tool_dir)
    if not current_version:
        return None

    dockerfile = tool_dir / current_version / "Dockerfile"
    if not dockerfile.exists():
        return None

    # --- 1. GitHub ---
    github_repo = extract_github_url(dockerfile, tool_name)
    if github_repo:
        upstream_tag = get_latest_release(github_repo, token)
        if not upstream_tag:
            upstream_tag = get_latest_tag(github_repo, token)
        if upstream_tag:
            return {
                "tool": tool_name,
                "current": current_version,
                "upstream": upstream_tag,
                "source": "github",
                "source_ref": github_repo,
                "outdated": is_newer(upstream_tag, current_version),
                "error": None,
            }
        return {
            "tool": tool_name,
            "current": current_version,
            "upstream": None,
            "source": "github",
            "source_ref": github_repo,
            "outdated": False,
            "error": "Could not fetch upstream version from GitHub",
        }

    # --- 2. GitLab ---
    gitlab_repo = extract_gitlab_url(dockerfile, tool_name)
    if gitlab_repo:
        upstream_tag = get_gitlab_latest_version(gitlab_repo)
        if upstream_tag:
            return {
                "tool": tool_name,
                "current": current_version,
                "upstream": upstream_tag,
                "source": "gitlab",
                "source_ref": gitlab_repo,
                "outdated": is_newer(upstream_tag, current_version),
                "error": None,
            }
        return {
            "tool": tool_name,
            "current": current_version,
            "upstream": None,
            "source": "gitlab",
            "source_ref": gitlab_repo,
            "outdated": False,
            "error": "Could not fetch upstream version from GitLab",
        }

    # --- 3. PyPI ---
    pypi_pkg = extract_pypi_package(dockerfile, tool_name)
    if pypi_pkg:
        upstream_version = get_pypi_latest_version(pypi_pkg)
        if upstream_version:
            return {
                "tool": tool_name,
                "current": current_version,
                "upstream": upstream_version,
                "source": "pypi",
                "source_ref": pypi_pkg,
                "outdated": is_newer(upstream_version, current_version),
                "error": None,
            }
        return {
            "tool": tool_name,
            "current": current_version,
            "upstream": None,
            "source": "pypi",
            "source_ref": pypi_pkg,
            "outdated": False,
            "error": "Could not fetch upstream version from PyPI",
        }

    # --- 4. Conda / Bioconda ---
    conda_info = extract_conda_package(dockerfile, tool_name)
    if conda_info:
        pkg_name, channel = conda_info
        upstream_version = get_conda_latest_version(pkg_name, channel)
        if upstream_version:
            return {
                "tool": tool_name,
                "current": current_version,
                "upstream": upstream_version,
                "source": "conda",
                "source_ref": f"{channel}/{pkg_name}",
                "outdated": is_newer(upstream_version, current_version),
                "error": None,
            }
        return {
            "tool": tool_name,
            "current": current_version,
            "upstream": None,
            "source": "conda",
            "source_ref": f"{channel}/{pkg_name}",
            "outdated": False,
            "error": f"Could not fetch upstream version from conda ({channel})",
        }

    # No upstream source found
    return {
        "tool": tool_name,
        "current": current_version,
        "upstream": None,
        "source": None,
        "source_ref": None,
        "outdated": False,
        "error": "No upstream source found in Dockerfile",
    }


def _source_link(result):
    """Return a Markdown link for the upstream source of a result dict."""
    source = result.get("source")
    ref = result.get("source_ref") or ""
    if source == "github":
        return f"[{ref}](https://github.com/{ref}) (GitHub)"
    if source == "gitlab":
        return f"[{ref}](https://gitlab.com/{ref}) (GitLab)"
    if source == "pypi":
        return f"[{ref}](https://pypi.org/project/{ref}/) (PyPI)"
    if source == "conda":
        channel, pkg = ref.split("/", 1) if "/" in ref else ("", ref)
        return f"[{ref}](https://anaconda.org/{channel}/{pkg}) (Conda)"
    return "N/A"


def main():
    token = os.environ.get("GITHUB_TOKEN")
    # Use GITHUB_REPOSITORY env var (set in GitHub Actions) for absolute URLs
    gh_repo_slug = os.environ.get("GITHUB_REPOSITORY", "evolbioinfo/dockerfiles")
    repo_base_url = f"https://github.com/{gh_repo_slug}/tree/main"

    tools = get_tool_dirs()

    results = []
    errors = []

    for tool in tools:
        result = check_tool(tool, token)
        if result is None:
            continue
        if result["error"]:
            errors.append(result)
        elif result["outdated"]:
            results.append(result)

    # Print Markdown summary
    if results:
        print("## Tools that could be updated\n")
        print("| Tool | Current Version | Latest Upstream | Source |")
        print("|------|----------------|-----------------|--------|")
        for r in sorted(results, key=lambda x: x["tool"]):
            tool_link = f"[{r['tool']}]({repo_base_url}/{r['tool']})"
            print(
                f"| {tool_link} "
                f"| `{r['current']}` "
                f"| `{r['upstream']}` "
                f"| {_source_link(r)} |"
            )
        print()
    else:
        print("## Tools that could be updated\n")
        print("All tools are up to date! 🎉\n")

    if errors:
        print("## Tools where the version could not be determined\n")
        print("| Tool | Current Version | Reason |")
        print("|------|----------------|--------|")
        for r in sorted(errors, key=lambda x: x["tool"]):
            tool_link = f"[{r['tool']}]({repo_base_url}/{r['tool']})"
            print(
                f"| {tool_link} "
                f"| `{r['current']}` "
                f"| {r['error']} |"
            )
        print()

    # Exit with non-zero status if there are outdated tools (useful in CI)
    if results:
        sys.exit(1)


if __name__ == "__main__":
    main()
