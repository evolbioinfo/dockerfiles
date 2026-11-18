# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repo is

A collection of Dockerfiles for bioinformatics tools (evolutionary biology / phylogenetics focus). Every image is published on Docker Hub as `evolbioinfo/<tool>:<version>`. There is no application code; the only scripts live in `.github/scripts/`.

## Layout

- `<tool>/<version>/Dockerfile`: one directory per published image tag. The directory name **is** the Docker tag (e.g. `raxml-ng/v2.0.2` → `evolbioinfo/raxml-ng:v2.0.2`). Exception: `bwa-mem/` is published as `evolbioinfo/bwa-mem2`.
- Version directory names are usually `v<semver>`. Unreleased builds use a commit hash: `dev<shorthash>` (e.g. `goalign/dev6b7572`), `<shorthash>`, or `<branch>-<hash>`. `dev/` and `latest/` are floating directories that track a branch.
- Some version directories hold extra files that the Dockerfile `COPY`s in (patches, wrapper `.sh` scripts). Files ending in `~` are editor backups and are gitignored.
- Base/environment images (`ubuntu`, `python`, `python-evol`, `r-base`, `r-evol`, …) are used by other images through `FROM evolbioinfo/<base>:<version>`. If you change one of them, the images built on it may need rebuilding.
- `ngphylogeny_multitools/` combines several tools in one image (e.g. `fasttree_gotree_goalign_parallel`).

## Adding a new image version

The usual workflow (see git history, e.g. "Added goalign dev6b7572"):

1. Copy the newest existing version directory of the tool to `<tool>/<newversion>/`.
2. Update the `# version` header comment and `ENV VERSION=...` (for commit builds, use the full commit hash in `VERSION`; `git checkout ${PREFIXVERSION}${VERSION}` is the common pattern). Check that the toolchain version (e.g. `GOVERSION` against upstream `go.mod`) is still compatible.
3. Build, smoke-test, push:
   ```bash
   cd <tool>/<version>
   docker build --platform linux/amd64 -t evolbioinfo/<tool>:<version> .
   docker run --rm --platform linux/amd64 evolbioinfo/<tool>:<version> --help   # or `version`
   docker push evolbioinfo/<tool>:<version>
   ```
   The images are amd64: many Dockerfiles download `linux-amd64` binaries, so pass `--platform linux/amd64` when building on Apple Silicon. Under Docker Desktop's Rosetta emulation, GNU tar in `ubuntu:26.04` fails with "Cannot open: Function not implemented": extract archives in an `ubuntu:24.04` stage and `COPY` the result (see `trimgalore/v2.5.0`), or build with Rosetta disabled (QEMU).
   For several images at once, `./build_updates.sh [--no-push] [--dry-run] [tool ...]` builds (and pushes) the `tool/version` entries of its `IMAGES` list, with one log per image in `build_logs/` (gitignored).
4. For new releases (not dev builds), update the "Latest Version" column for the tool in the `README.md` tables.
5. Commit only the new directory with `git add <tool>/<version>`. The working tree usually has many untracked scratch files (logs, test outputs, `UPDATE*.md`, `run_update.sh`), so do not use `git add -A`.

## Dockerfile conventions (from README "Contributing")

- `LABEL maintainer="..."` (older files use the deprecated `MAINTAINER`).
- `ENV VERSION=<version>` holds the tool version.
- Install, build, and clean up in a single `RUN` step: remove build deps, `apt-get autoremove`, `apt-get clean`, `rm -rf /var/lib/apt/lists/*`, and delete toolchains/sources (e.g. `/gopath /usr/local/go`).
- End the build with `mkdir /pasteur` (the shared mount point convention).
- Set `ENTRYPOINT` to the tool's executable where it makes sense.
- Go tools (goalign, gotree) build with `CGO_ENABLED=0` and run `make build && make test && make testcommands` during the image build, so a successful build also means the upstream tests passed.

## Automation (`.github/`)

- `check_updates.py` (weekly workflow `check_updates.yml`) picks the latest version directory of each tool (it ignores `dev`/`latest` and sorts semver directories above commit hashes), reads the Dockerfile to work out the upstream source (GitHub/GitLab/PyPI/Conda/SourceForge), and posts the outdated tools to a GitHub issue labelled `update-check`. It only considers stable upstream versions (it ignores pre-releases). To run it locally: `GITHUB_TOKEN=$(gh auth token) python3 .github/scripts/check_updates.py > UPDATE.md` (needs `packaging`). The report goes to stdout and warnings to stderr; the exit code is 0 unless the script crashes. Run locally, it also sees untracked version directories, so its "current version" can differ from the CI run. Tools that use commit or custom versioning are listed in `SKIP_VERSION_CHECK`; add new ones there. For a tool whose Dockerfile doesn't mention a supported source (e.g. it downloads from the tool's own website), add an entry to `UPSTREAM_OVERRIDES`. For a GitHub project without releases or tags, a `github_file` override reads the version from a source file (`path` + a regex `pattern` whose first group is the version).
- `BIOTOOLS.md` (the tool → bio.tools ID mapping) is a static file and nothing regenerates it any more.
