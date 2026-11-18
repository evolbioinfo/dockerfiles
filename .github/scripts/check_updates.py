#!/usr/bin/env python3
"""
Check for newer versions of tools in the evolbioinfo/dockerfiles repository.

For each tool, this script:
1. Finds the latest version directory in the tool's folder
2. Determines the upstream source: UPSTREAM_OVERRIDES first, otherwise it is
   detected from the Dockerfile (GitHub, GitLab, PyPI, Conda, SourceForge).
   Overrides can also read the version from a file of a GitHub repository
   ('github_file'), for projects without releases or tags
3. Queries the upstream for its latest *stable* version (pre-releases ignored)
4. Compares it to the current version and reports outdated tools

The Markdown report is written to stdout and warnings to stderr.
Exit status: 0 when the check ran (whether or not tools are outdated),
non-zero if the script itself crashed.
"""

import json
import os
import re
import sys
import traceback
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from packaging.version import InvalidVersion, Version


REPO_ROOT = Path(__file__).parent.parent.parent

# Directories to skip (not real tools)
SKIP_DIRS = {".git", ".github"}

# Directories that use non-release versioning (commits, custom schemes)
SKIP_VERSION_CHECK = {
    "artic-ncov2019",  # commit hash
    "ml_bootstrap",  # commit hash
    "ngphylogeny_multitools",  # branch name
    "pcoc",  # commit hash
    "polecat",  # commit hash
    "ptp",  # commit hash
    "reseq",  # commit hash
    "snag",  # 'master'
    "strainline",  # commit hash
    "table2itol",  # commit hash
    "treedater",  # commit hash
    "treestructure",  # commit hash
    "inkscape",  # 'latest' tag
    "sdrmhunter",  # custom script; pip deps are not the tool itself
    # No upstream source that can be queried (no repository, package or API)
    "guppy",  # Oxford Nanopore, closed source and discontinued (replaced by dorado)
    "jphmm",  # tarball on jphmm.gobics.de, custom version (03.2015)
    "nextstrain-base",  # Docker Hub build tags (build-<timestamp>)
    "noisy",  # tarball on bioinf.uni-leipzig.de
    "vivan",  # zip on vivan.tau.ac.il
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

# Upstream sources for tools whose Dockerfile does not reference a supported
# source (e.g. the tool is downloaded from its own website). These take
# precedence over detection from the Dockerfile.
UPSTREAM_OVERRIDES = {
    "admixture": {"source": "conda", "ref": "admixture", "channels": ["bioconda"]},
    # No releases nor tags: the version is read from a source file
    # ('github_file': first regex group of 'pattern' in 'path', default branch)
    "bayestraits": {
        "source": "github_file", "ref": "AndrewPMeade/BayesTraits-Release",
        "path": "src/TypeDef.h", "pattern": r"BayesTraits V(\d+(?:\.\d+)+)",
    },
    "newick_utilities": {
        "source": "github_file", "ref": "tjunier/newick_utils",
        "path": "configure.ac", "pattern": r"AC_INIT\(\[[^]]*\],\s*\[([^]]+)\]",
    },
    "papara": {
        "source": "github_file", "ref": "sim82/papara_nt",
        "path": "papara.cpp", "pattern": r'get_version_string\(\)\s*\{\s*return\s+std::string\(\s*"([^"]+)"',
    },
    "clustal_omega": {"source": "conda", "ref": "clustalo", "channels": ["bioconda"]},
    "deseq": {"source": "conda", "ref": "bioconductor-deseq", "channels": ["bioconda"]},
    "finestructure": {"source": "conda", "ref": "finestructure", "channels": ["bioconda"]},
    "gblocks": {"source": "conda", "ref": "gblocks", "channels": ["bioconda"]},
    "hmmer": {"source": "github", "ref": "EddyRivasLab/hmmer"},
    "igv": {"source": "github", "ref": "igvteam/igv"},
    "label": {"source": "github", "ref": "CDCgov/label"},
    "lastal": {"source": "gitlab", "host": "gitlab.com", "ref": "mcfrith/last"},
    "multiqc": {"source": "pypi", "ref": "multiqc"},
    "paml": {"source": "github", "ref": "abacus-gene/paml"},
    "phyml-sms": {"source": "gitlab", "host": "gite.lirmm.fr", "ref": "atgc/sms"},
    "predicthaplo": {"source": "github", "ref": "cbg-ethz/PredictHaplo"},
    "sratoolkit": {"source": "github", "ref": "ncbi/sra-tools"},
    "stringtie": {"source": "github", "ref": "gpertea/stringtie"},
    "tqdist": {"source": "conda", "ref": "tqdist", "channels": ["bioconda"]},
    "velvet": {"source": "github", "ref": "dzerbino/velvet"},
    "wget": {"source": "conda", "ref": "wget", "channels": ["conda-forge"]},
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

# Conda channels tried, in order, after any channel given explicitly in the
# Dockerfile. Names are anaconda.org owners.
CONDA_CHANNEL_PRIORITY = ["bioconda", "conda-forge", "anaconda"]
# The 'defaults' channel is published on anaconda.org as 'anaconda'.
CONDA_CHANNEL_ALIASES = {"defaults": "anaconda", "main": "anaconda"}

SOURCE_LABELS = {
    "github": "GitHub",
    "gitlab": "GitLab",
    "pypi": "PyPI",
    "conda": "Conda",
    "sourceforge": "SourceForge",
    "github_file": "GitHub (version in source file)",
}

USER_AGENT = "evolbioinfo-dockerfiles-checker"
HTTP_TIMEOUT = 15
MAX_WORKERS = 8
GITHUB_MAX_TAG_PAGES = 5  # 100 tags per page


# ---------------------------------------------------------------------------
# Version parsing
# ---------------------------------------------------------------------------

# Git commit hashes, optionally prefixed with 'v' or 'dev' (e.g. v4bb2daf,
# dev6b7572, 581c23199). Must contain both a digit and a hex letter.
_COMMIT_RE = re.compile(
    r"(?:v|dev)?(?=[0-9a-f]*[a-f])(?=[0-9a-f]*[0-9])[0-9a-f]{6,40}", re.IGNORECASE
)
# A dotted version after an optional non-numeric prefix: v1.2, v.2.4.1,
# MAPLEv0.7.5, jq-1.8.2, Version_13.46.2.7c9e712d, v4.4-v1.41-v0.3.19
_DOTTED_RE = re.compile(r"^[A-Za-z_.-]*(\d+(?:\.\d+)+)(.*)$")
# A single number (v980, 1654, r1613), only with at most a v/r prefix and
# nothing numeric after it, so that tags such as 'test9' or 'VDB-4760-2.9.7'
# are not read as versions 9 or 4760.
_PLAIN_RE = re.compile(r"^(?:[vVrR][-_.]?)?(\d+)(\D*)$")
# Pre-release suffix right after the number: -alpha2, -RC2, b3, -beta.8, .dev
_PRE_RE = re.compile(
    r"^[-_.]?(alpha|beta|preview|pre|rc|dev|a|b|c)[-_.]?(\d*)(?![A-Za-z0-9])",
    re.IGNORECASE,
)
_PRE_TAGS = {
    "alpha": "a", "a": "a",
    "beta": "b", "b": "b",
    "rc": "rc", "c": "rc", "pre": "rc", "preview": "rc",
    "dev": ".dev",
}


def is_commit_hash(name):
    return bool(_COMMIT_RE.fullmatch(name))


def parse_version(text):
    """
    Parse a version from a directory name or upstream tag.

    Returns a packaging Version, or None if no version can be found.
    A single trailing letter (v0.3.6a, 4.8a) is a revision of that version and
    is ordered after it; alpha/beta/rc suffixes are pre-releases.
    """
    if not text or is_commit_hash(text):
        return None
    m = _DOTTED_RE.match(text) or _PLAIN_RE.match(text)
    if not m:
        return None
    number, rest = m.group(1), m.group(2)
    suffix = ""
    if re.fullmatch(r"[A-Za-z]", rest):
        suffix = f".post{ord(rest.lower()) - ord('a') + 1}"
    else:
        pre = _PRE_RE.match(rest)
        if pre:
            suffix = f"{_PRE_TAGS[pre.group(1).lower()]}{pre.group(2) or '0'}"
    try:
        return Version(number + suffix)
    except InvalidVersion:
        return None


def is_stable(text):
    v = parse_version(text)
    return v is not None and not v.is_prerelease


def latest_stable(names):
    """Return the name with the highest stable version, or None."""
    best, best_v = None, None
    for name in names:
        if not name:
            continue
        v = parse_version(name)
        if v is None or v.is_prerelease:
            continue
        if best_v is None or v > best_v:
            best, best_v = name, v
    return best


def _version_dir_sort_key(name):
    """
    Sort key for version directories: parseable versions sort above commit
    hashes and unparseable names. Ties (e.g. pangolin's
    v4.4-v1.38-... vs v4.4-v1.41-...) are broken on all numbers in the name.
    """
    v = parse_version(name)
    if v is None:
        return (0, Version("0"), (), name)
    numbers = tuple(int(n) for n in re.findall(r"\d+", name))
    return (1, v, numbers, name)


def get_latest_version_dir(tool_dir):
    """
    Return the latest versioned subdirectory for a tool.
    Skips 'dev' and 'latest' directories.
    """
    version_dirs = [
        e.name for e in tool_dir.iterdir()
        if e.is_dir() and e.name not in ("dev", "latest")
    ]
    if not version_dirs:
        return None
    return max(version_dirs, key=_version_dir_sort_key)


def is_newer(upstream, current):
    """
    Return True if upstream is newer than current, False if not, and None if
    either version cannot be parsed.
    """
    upstream_v = parse_version(upstream)
    current_v = parse_version(current)
    if upstream_v is None or current_v is None:
        return None
    return upstream_v > current_v


# ---------------------------------------------------------------------------
# Dockerfile parsing
# ---------------------------------------------------------------------------

def _name_key(name):
    return name.lower().replace("-", "").replace("_", "")


def _matches_tool(name, tool_name):
    if not tool_name:
        return False
    a, b = _name_key(name), _name_key(tool_name)
    return a in b or b in a


def _pick_for_tool(candidates, tool_name, name_of=lambda c: c):
    """Return the first candidate whose name matches the tool, or None."""
    for candidate in candidates:
        if _matches_tool(name_of(candidate), tool_name):
            return candidate
    return None


def _comment_lines(content):
    return [line for line in content.splitlines() if line.startswith("#")]


def _shell_commands(content):
    """
    Yield the individual shell commands of a Dockerfile, with line
    continuations joined so that multi-line installs are seen as one command.
    """
    joined = re.sub(r"\\[ \t]*\r?\n", " ", content)
    for line in joined.splitlines():
        if line.lstrip().startswith("#"):
            continue
        for command in re.split(r"&&|\|\||;|\|", line):
            yield command


_GITHUB_URL_RE = re.compile(r"https://github\.com/([A-Za-z0-9_.-]+)/([A-Za-z0-9_.-]+)")


def extract_github_repo(content, tool_name=None):
    """
    Return 'owner/repo' for the tool's GitHub repository, or None.

    Header comments are trusted (preferring a repo named like the tool);
    URLs in the body are only used if the repo name matches the tool, to
    avoid picking up bundled dependencies.
    """
    def repos(lines):
        found = []
        for line in lines:
            for m in _GITHUB_URL_RE.finditer(line):
                found.append(f"{m.group(1)}/{m.group(2).removesuffix('.git')}")
        return found

    def repo_name(repo):
        return repo.split("/")[1]

    header = repos(_comment_lines(content))
    if header:
        return _pick_for_tool(header, tool_name, repo_name) or header[0]
    return _pick_for_tool(repos(content.splitlines()), tool_name, repo_name)


_GITLAB_URL_RE = re.compile(r"https?://((?:gitlab|gite)[A-Za-z0-9.-]*)/([A-Za-z0-9_.~/-]+)")
# Path segments that end the project path in a GitLab URL
_GITLAB_PATH_STOP = {
    "-", "raw", "blob", "tree", "repository", "uploads", "archive",
    "tags", "releases", "commits", "commit", "wikis",
}


def _clean_gitlab_path(raw):
    """Return the 'group[/subgroup]/project' part of a GitLab URL path."""
    parts = []
    for segment in raw.split("/"):
        if not segment or segment in _GITLAB_PATH_STOP:
            break
        parts.append(segment)
    if len(parts) < 2:
        return None
    parts[-1] = parts[-1].removesuffix(".git")
    return "/".join(parts)


def extract_gitlab_project(content, tool_name=None):
    """
    Return (host, 'namespace/project') for a gitlab.com or self-hosted GitLab
    project, or None. Same trust rules as extract_github_repo.
    """
    def projects(lines):
        found = []
        for line in lines:
            for m in _GITLAB_URL_RE.finditer(line):
                path = _clean_gitlab_path(m.group(2))
                if path:
                    found.append((m.group(1), path))
        return found

    def project_name(candidate):
        return candidate[1].rsplit("/", 1)[-1]

    header = projects(_comment_lines(content))
    if header:
        return _pick_for_tool(header, tool_name, project_name) or header[0]
    return _pick_for_tool(projects(content.splitlines()), tool_name, project_name)


_PIP_INSTALL_RE = re.compile(r"\bpip(?:3(?:\.\d+)?)?\s+install\b(.*)")
_PIP_PKG_RE = re.compile(
    r"(?<![\w.-])([A-Za-z0-9][A-Za-z0-9_.-]*)(?:\[[^\]]*\])?=="
    r"(\$\{?[A-Za-z_]\w*\}?|\d[\w.+!-]*)"
)


def extract_pypi_package(content, tool_name=None):
    """
    Return the PyPI package from a `pip install <pkg>==<version>` command
    that best matches tool_name, or None.
    """
    packages = []
    for command in _shell_commands(content):
        m = _PIP_INSTALL_RE.search(command)
        if not m:
            continue
        for pm in _PIP_PKG_RE.finditer(m.group(1)):
            if pm.group(1).lower() not in COMMON_UTILITY_PACKAGES:
                packages.append(pm.group(1))
    if not packages:
        return None
    return _pick_for_tool(packages, tool_name) or packages[0]


_CONDA_INSTALL_RE = re.compile(
    r"\b(?:micromamba|mamba|conda)\s+(?:install|create)\b(.*)", re.IGNORECASE
)
_CONDA_CHANNEL_RE = re.compile(r"(?:^|\s)(?:-c|--channel)(?:\s+|=)(\S+)")
_CONDA_PKG_RE = re.compile(
    r"(?<![\w.:-])(?:([A-Za-z0-9][\w.-]*)::)?([A-Za-z0-9][\w.-]*)==?"
    r"(\$\{?[A-Za-z_]\w*\}?|\d[\w.*+-]*)"
)


def _normalize_channel(channel):
    name = channel.strip("'\"").rstrip("/").rsplit("/", 1)[-1]
    return CONDA_CHANNEL_ALIASES.get(name, name)


def _channel_order(channels):
    """Explicit channels first (in priority order), then the defaults."""
    def rank(ch):
        return CONDA_CHANNEL_PRIORITY.index(ch) if ch in CONDA_CHANNEL_PRIORITY else len(CONDA_CHANNEL_PRIORITY)
    ordered = sorted(dict.fromkeys(channels), key=rank)
    return ordered + [ch for ch in CONDA_CHANNEL_PRIORITY if ch not in ordered]


def extract_conda_package(content, tool_name=None):
    """
    Return (package, [channels to try]) from a conda/mamba install command,
    or None.
    """
    candidates = []
    for command in _shell_commands(content):
        m = _CONDA_INSTALL_RE.search(command)
        if not m:
            continue
        args = m.group(1)
        channels = [_normalize_channel(c) for c in _CONDA_CHANNEL_RE.findall(args)]
        for pm in _CONDA_PKG_RE.finditer(args):
            pkg = pm.group(2)
            if pkg.lower() in COMMON_UTILITY_PACKAGES:
                continue
            pkg_channels = [_normalize_channel(pm.group(1))] if pm.group(1) else []
            candidates.append((pkg, list(dict.fromkeys(pkg_channels + _channel_order(channels)))))
    if not candidates:
        return None
    return _pick_for_tool(candidates, tool_name, lambda c: c[0]) or candidates[0]


_SOURCEFORGE_RES = [
    re.compile(r"sourceforge\.net/projects/([A-Za-z0-9_.-]+)"),
    re.compile(r"downloads\.sourceforge\.net/project/([A-Za-z0-9_.-]+)"),
    re.compile(r"svn\.code\.sf\.net/p/([A-Za-z0-9_.-]+)"),
    re.compile(r"(?:^|[/\s(\"'])([A-Za-z0-9_.-]+)\.sourceforge\.net"),
]


def extract_sourceforge_project(content, tool_name=None):
    """
    Return the SourceForge project name, or None.

    Only returns a project whose name matches the tool, to avoid picking up
    SourceForge URLs of bundled dependencies.
    """
    header, body = [], []
    for line in content.splitlines():
        target = header if line.strip().startswith("#") else body
        for pattern in _SOURCEFORGE_RES:
            for m in pattern.finditer(line):
                target.append(m.group(1).rstrip("/"))
    candidates = header + body
    if not candidates:
        return None
    if tool_name:
        return _pick_for_tool(candidates, tool_name)
    return candidates[0]


def detect_source(tool_name, content):
    """Return a source dict ({'source': ..., 'ref': ..., ...}) or None."""
    if tool_name in UPSTREAM_OVERRIDES:
        return dict(UPSTREAM_OVERRIDES[tool_name])
    repo = extract_github_repo(content, tool_name)
    if repo:
        return {"source": "github", "ref": repo}
    gitlab = extract_gitlab_project(content, tool_name)
    if gitlab:
        return {"source": "gitlab", "host": gitlab[0], "ref": gitlab[1]}
    pkg = extract_pypi_package(content, tool_name)
    if pkg:
        return {"source": "pypi", "ref": pkg}
    conda = extract_conda_package(content, tool_name)
    if conda:
        return {"source": "conda", "ref": conda[0], "channels": conda[1]}
    project = extract_sourceforge_project(content, tool_name)
    if project:
        return {"source": "sourceforge", "ref": project}
    return None


# ---------------------------------------------------------------------------
# Upstream queries
# ---------------------------------------------------------------------------

class UpstreamError(Exception):
    """The upstream could not be queried (network error, HTTP error, ...)."""


def http_get(url, headers=None):
    """Return the response body as text, None on HTTP 404, or raise UpstreamError."""
    req = urllib.request.Request(url)
    req.add_header("User-Agent", USER_AGENT)
    for key, value in (headers or {}).items():
        req.add_header(key, value)
    try:
        with urllib.request.urlopen(req, timeout=HTTP_TIMEOUT) as resp:
            return resp.read().decode(errors="replace")
    except urllib.error.HTTPError as e:
        if e.code == 404:
            return None
        hint = " (rate limited? set GITHUB_TOKEN)" if e.code in (403, 429) and "api.github.com" in url else ""
        raise UpstreamError(f"HTTP {e.code}{hint}") from e
    except (urllib.error.URLError, OSError) as e:
        raise UpstreamError(f"network error: {e}") from e


def http_get_json(url, headers=None):
    body = http_get(url, headers)
    if body is None:
        return None
    try:
        return json.loads(body)
    except json.JSONDecodeError as e:
        raise UpstreamError(f"invalid JSON from {url}") from e


def _github_headers(token):
    headers = {
        "Accept": "application/vnd.github+json",
        "X-GitHub-Api-Version": "2022-11-28",
    }
    if token:
        headers["Authorization"] = f"Bearer {token}"
    return headers


def github_latest(repo, token=None):
    """
    Latest stable version of a GitHub repo: the latest release if it is a
    stable version, otherwise the highest stable tag. The tags API is not
    sorted by version, so several pages are scanned.
    """
    headers = _github_headers(token)
    release = http_get_json(f"https://api.github.com/repos/{repo}/releases/latest", headers)
    if release and not release.get("prerelease") and is_stable(release.get("tag_name")):
        return release["tag_name"]
    names = []
    for page in range(1, GITHUB_MAX_TAG_PAGES + 1):
        tags = http_get_json(
            f"https://api.github.com/repos/{repo}/tags?per_page=100&page={page}", headers
        )
        if not tags:
            break
        names.extend(t["name"] for t in tags)
        if len(tags) < 100:
            break
    return latest_stable(names)


def gitlab_latest(host, project):
    """Latest stable version of a GitLab project (releases, then tags)."""
    base = f"https://{host}/api/v4/projects/{urllib.parse.quote(project, safe='')}"
    releases = http_get_json(f"{base}/releases?per_page=20") or []
    tag = latest_stable(
        r.get("tag_name") for r in releases if not r.get("upcoming_release")
    )
    if tag:
        return tag
    tags = http_get_json(f"{base}/repository/tags?per_page=100") or []
    return latest_stable(t.get("name") for t in tags)


def pypi_latest(package):
    """Latest stable version of a PyPI package."""
    data = http_get_json(f"https://pypi.org/pypi/{package}/json")
    if not data:
        return None
    releases = data.get("releases") or {}
    available = [v for v, files in releases.items() if not files or not all(f.get("yanked") for f in files)]
    return latest_stable(available) or latest_stable([data.get("info", {}).get("version")])


def conda_latest(package, channels):
    """Return (latest stable version, channel) from the first channel that has the package."""
    for channel in channels:
        data = http_get_json(f"https://api.anaconda.org/package/{channel}/{package}")
        if data:
            version = latest_stable(data.get("versions") or []) or latest_stable([data.get("latest_version")])
            return version, channel
    return None, None


def sourceforge_latest(project):
    """
    Latest stable version of a SourceForge project, from the version numbers
    embedded in the file names of its release RSS feed.
    """
    raw = http_get(f"https://sourceforge.net/projects/{project}/rss?limit=100")
    if raw is None:
        return None
    try:
        root = ET.fromstring(raw)
    except ET.ParseError as e:
        raise UpstreamError(f"invalid RSS feed: {e}") from e

    # BBMap_39.06.tar.gz → 39.06, musket-1.1.tar.gz → 1.1,
    # subread-2.0.4-source.tar.gz → 2.0.4
    version_re = re.compile(r"[-_v](\d+(?:\.\d+)+)", re.IGNORECASE)
    versions = []
    for item in root.iter("item"):
        title = item.find("title")
        if title is None or not title.text:
            continue
        # RSS titles are paths like "/subdir/BBMap_39.06.tar.gz"; take basename
        filename = title.text.strip().rstrip("/").rsplit("/", 1)[-1]
        if "." not in filename:
            continue  # directory
        m = version_re.search(filename)
        if m:
            versions.append(m.group(1))
    return latest_stable(versions)


def github_file_version(repo, path, pattern):
    """
    Version read from a file on the default branch of a GitHub repository, for
    projects that publish neither releases nor tags. 'pattern' is a regex
    whose first group is the version.
    """
    content = http_get(f"https://raw.githubusercontent.com/{repo}/HEAD/{path}")
    if content is None:
        raise UpstreamError(f"{path} not found in {repo}")
    m = re.search(pattern, content)
    if not m:
        raise UpstreamError(f"no version matching {pattern!r} in {path}")
    return m.group(1)


def fetch_latest(src, token=None):
    """Return the latest stable upstream version for a source dict, or None."""
    kind = src["source"]
    if kind == "github":
        return github_latest(src["ref"], token)
    if kind == "gitlab":
        return gitlab_latest(src["host"], src["ref"])
    if kind == "pypi":
        return pypi_latest(src["ref"])
    if kind == "conda":
        version, channel = conda_latest(src["ref"], src["channels"])
        src["channel"] = channel or src["channels"][0]
        return version
    if kind == "sourceforge":
        return sourceforge_latest(src["ref"])
    if kind == "github_file":
        return github_file_version(src["ref"], src["path"], src["pattern"])
    raise ValueError(f"unknown source {kind!r}")


def source_ref(src):
    if src["source"] == "conda":
        return f"{src.get('channel') or src['channels'][0]}/{src['ref']}"
    return src["ref"]


def source_url(src):
    kind, ref = src["source"], src["ref"]
    if kind == "github":
        return f"https://github.com/{ref}"
    if kind == "gitlab":
        return f"https://{src['host']}/{ref}"
    if kind == "pypi":
        return f"https://pypi.org/project/{ref}/"
    if kind == "conda":
        return f"https://anaconda.org/{source_ref(src)}"
    if kind == "sourceforge":
        return f"https://sourceforge.net/projects/{ref}/"
    if kind == "github_file":
        return f"https://github.com/{ref}/blob/HEAD/{src['path']}"
    return None


# ---------------------------------------------------------------------------
# Checking and reporting
# ---------------------------------------------------------------------------

def check_tool(tool_name, token=None):
    """
    Check a single tool for available updates.

    Returns None if the tool is not checked, otherwise a dict with:
      - tool, current: tool name and current version directory
      - upstream: latest stable upstream version (or None)
      - source, source_ref, source_url: where the upstream version comes from
      - outdated: True if upstream > current
      - error: why the tool could not be checked (or None)
    """
    if tool_name in SKIP_VERSION_CHECK:
        return None
    tool_dir = REPO_ROOT / tool_name
    if not tool_dir.is_dir():
        return None
    current = get_latest_version_dir(tool_dir)
    if not current:
        return None
    dockerfile = tool_dir / current / "Dockerfile"
    try:
        content = dockerfile.read_text(errors="replace")
    except OSError:
        return None

    result = {
        "tool": tool_name,
        "current": current,
        "upstream": None,
        "source": None,
        "source_ref": None,
        "source_url": None,
        "outdated": False,
        "error": None,
    }

    src = detect_source(tool_name, content)
    if src is None:
        result["error"] = "No upstream source found in Dockerfile (see UPSTREAM_OVERRIDES)"
        return result

    label = SOURCE_LABELS[src["source"]]
    result["source"] = src["source"]
    try:
        upstream = fetch_latest(src, token)
    except UpstreamError as e:
        upstream = None
        result["error"] = f"Could not query {label}: {e}"
    result["source_ref"] = source_ref(src)
    result["source_url"] = source_url(src)
    if result["error"]:
        return result
    if not upstream:
        result["error"] = f"No stable release or tag found on {label}"
        return result

    result["upstream"] = upstream
    newer = is_newer(upstream, current)
    if newer is None:
        result["error"] = f"Cannot compare versions `{current}` and `{upstream}`"
    else:
        result["outdated"] = newer
    return result


def safe_check_tool(tool_name, token=None):
    """check_tool, turning unexpected exceptions into an error row."""
    try:
        return check_tool(tool_name, token)
    except Exception as e:  # keep checking the other tools
        traceback.print_exc(file=sys.stderr)
        return {
            "tool": tool_name, "current": "?", "upstream": None,
            "source": None, "source_ref": None, "source_url": None,
            "outdated": False, "error": f"Unexpected error: {type(e).__name__}: {e}",
        }


def get_tool_dirs():
    """Return a sorted list of tool directory names."""
    return sorted(
        e.name for e in REPO_ROOT.iterdir()
        if e.is_dir() and e.name not in SKIP_DIRS and not e.name.startswith(".")
    )


def _md_cell(text):
    return str(text).replace("|", "\\|")


def _source_link(result):
    """Return a Markdown link for the upstream source of a result dict."""
    if not result.get("source"):
        return "N/A"
    label = SOURCE_LABELS[result["source"]]
    return f"[{result['source_ref']}]({result['source_url']}) ({label})"


def main():
    token = os.environ.get("GITHUB_TOKEN")
    if not token:
        print("Warning: GITHUB_TOKEN is not set; GitHub API calls may be rate limited.", file=sys.stderr)
    # Set by GitHub Actions; defaults for local runs
    repo_slug = os.environ.get("GITHUB_REPOSITORY", "evolbioinfo/dockerfiles")
    branch = os.environ.get("GITHUB_REF_NAME") or "master"
    repo_base_url = f"https://github.com/{repo_slug}/tree/{branch}"

    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as pool:
        checked = list(pool.map(lambda t: safe_check_tool(t, token), get_tool_dirs()))
    checked = [r for r in checked if r is not None]
    outdated = [r for r in checked if r["outdated"]]
    errors = [r for r in checked if r["error"]]

    def tool_link(r):
        return f"[{r['tool']}]({repo_base_url}/{r['tool']})"

    print("## Tools that could be updated\n")
    if outdated:
        print("| Tool | Current Version | Latest Upstream | Source |")
        print("|------|----------------|-----------------|--------|")
        for r in outdated:
            print(f"| {tool_link(r)} | `{r['current']}` | `{r['upstream']}` | {_source_link(r)} |")
    else:
        print("All tools are up to date! 🎉")
    print()

    if errors:
        print("## Tools where the version could not be determined\n")
        print("| Tool | Current Version | Source | Reason |")
        print("|------|----------------|--------|--------|")
        for r in errors:
            print(f"| {tool_link(r)} | `{r['current']}` | {_source_link(r)} | {_md_cell(r['error'])} |")
        print()


if __name__ == "__main__":
    main()
