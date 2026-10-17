"""Check that every extra this package and its docs name is one pip installs.

pip treats an extra that a package does not declare as a warning, not an
error. Up to 0.15.0 this package declared none, so ``pip install "yamc[cad]"``
printed ``yamc 0.15.0 does not provide the extra 'cad'``, installed yamc-core
without cadquery and exited 0, while installation.md, geometry.md and
yamc-core's own error messages kept telling people to run it. The extras
existed only on yamc-core.

So this checks the three places an extra has to agree:

* pyproject.toml forwards each extra to the yamc-core extra of the same name,
  at the same exact pin as the plain dependency. An extra pinned apart from it
  could not be resolved.
* The pinned yamc-core release, as PyPI publishes it, declares every extra
  forwarded here, and every extra it declares is forwarded, except those in
  NOT_FORWARDED.
* Every ``yamc[...]`` the README and the docs pages name is declared here, and
  every ``yamc-core[...]`` they name is declared by that yamc-core release.

It reads PyPI's JSON API rather than an installed wheel, so it needs no install
and sees the metadata pip will see. The publish workflow runs it as well, where
it is also the check that the pinned yamc-core is on PyPI at all.

Run: python docs/check_extras.py
"""

from __future__ import annotations

import json
import re
import sys
import tomllib
import urllib.error
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

# yamc-core extras this package does not forward, and why.
NOT_FORWARDED = {
    "dev": "the tooling for building yamc-core from a checkout of fusion-neutronics/core",
}

# The plain dependency and each extra's requirement: `yamc-core==0.16.0`,
# `yamc-core[cad]==0.16.0`.
PIN = re.compile(r"^yamc-core(?:\[([^\]]+)\])?==([^\s;]+)$")

# `yamc[cad]`, `yamc[all,pandas]` or `yamc-core[all]` anywhere in a page, in
# prose or in a command. An extra name starts with a letter, which keeps
# indexing such as `yamc[0]` out.
MENTION = re.compile(r"\b(yamc(?:-core)?)\[\s*([A-Za-z][\w.\-]*(?:\s*,\s*[A-Za-z][\w.\-]*)*)\s*\]")


def normalize(extra: str) -> str:
    """Compare extra names the way pip does (PEP 685)."""
    return re.sub(r"[-_.]+", "-", extra).lower()


def main() -> int:
    project = tomllib.loads((ROOT / "pyproject.toml").read_text())["project"]
    failures: list[str] = []

    deps = project.get("dependencies", [])
    pin = PIN.match(deps[0]) if len(deps) == 1 else None
    if not pin or pin.group(1):
        print(f"::error::pyproject.toml must depend on exactly one exact yamc-core pin, not {deps}")
        return 1
    version = pin.group(2)
    print(f"pinned yamc-core=={version}")

    extras = {normalize(name): reqs for name, reqs in project.get("optional-dependencies", {}).items()}
    for name, reqs in extras.items():
        want = f"yamc-core[{name}]=={version}"
        if reqs != [want]:
            failures.append(f"pyproject.toml: the extra '{name}' is {reqs}, not ['{want}']")

    url = f"https://pypi.org/pypi/yamc-core/{version}/json"
    try:
        with urllib.request.urlopen(url, timeout=60) as response:
            info = json.load(response)["info"]
    except urllib.error.HTTPError as e:
        print(f"::error::yamc-core {version} is not on PyPI (HTTP {e.code}). Release yamc-core first.")
        return 1
    if info.get("provides_extra") is None:
        # Not an empty list: metadata that says nothing. Guessing from
        # requires_dist would miss an extra with no requirements of its own.
        print(f"::error::PyPI lists no Provides-Extra for yamc-core {version}, so its extras cannot be checked")
        return 1
    provided = {normalize(name) for name in info["provides_extra"]}
    print(f"yamc-core {version} declares: {', '.join(sorted(provided))}")

    for name in sorted(extras.keys() - provided):
        failures.append(f"pyproject.toml: the extra '{name}' forwards to yamc-core[{name}], "
                        f"which yamc-core {version} does not declare")
    for name in sorted(provided - extras.keys() - NOT_FORWARDED.keys()):
        failures.append(f"pyproject.toml: yamc-core {version} declares the extra '{name}', "
                        f"which is not forwarded (add it, or list it in NOT_FORWARDED with the reason)")

    pages = [ROOT / "README.md", *sorted((ROOT / "docs").rglob("*.md"))]
    declared = {"yamc": (project["version"], extras.keys()), "yamc-core": (version, provided)}
    named = 0
    for page in pages:
        for number, line in enumerate(page.read_text().splitlines(), 1):
            for package, names in MENTION.findall(line):
                release, known = declared[package]
                for name in (normalize(n.strip()) for n in names.split(",")):
                    named += 1
                    if name not in known:
                        failures.append(f"{page.relative_to(ROOT)}:{number}: {package}[{name}] is not an extra "
                                        f"{package} {release} declares")
    print(f"{named} extras named across {len(pages)} pages")

    if failures:
        print("::error::extras that `pip install` would not honour:", *failures, sep="\n  ")
        return 1
    print(f"forwarded: {', '.join(sorted(extras))}; not forwarded: {', '.join(sorted(NOT_FORWARDED))}: OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
