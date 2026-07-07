#!/usr/bin/env python3
"""Verify that the GPU coverage audit tracks non-LVLSET lesgo.conf keys."""

from __future__ import annotations

import argparse
import re
import sys

from repo_paths import repo_path

SOURCE_PATH = repo_path("input_util.f90")
DOC_PATH = repo_path("docs", "gpu_port_coverage_audit.md")

BLOCK_SUBROUTINES = {
    "DOMAIN": "domain_block",
    "MODEL": "model_block",
    "CORIOLIS": "coriolis_block",
    "FLOW_COND": "flow_cond_block",
    "OUTPUT": "output_block",
    "TURBINES": "turbines_block",
    "SCALARS": "scalars_block",
}

HIT_KEYS = {
    "UP_IN",
    "TI_OUT",
    "LX_HIT",
    "LY_HIT",
    "LZ_HIT",
    "NX_HIT",
    "NY_HIT",
    "NZ_HIT",
    "U_FILE",
    "V_FILE",
    "W_FILE",
}

CASE_KEY_RE = re.compile(r"case\s*\(\s*'(?P<key>[A-Z0-9_]+)'\s*\)", re.IGNORECASE)
ROW_RE = re.compile(r"^\|\s*`(?P<block>[A-Z_]+)`\s*\|\s*(?P<keys>.*?)\s*\|")
KEY_RE = re.compile(r"`(?P<key>[A-Z0-9_]+)`")


def subroutine_text(source: str, name: str) -> str:
    pattern = re.compile(
        rf"subroutine\s+{re.escape(name)}\s*\(\)(?P<body>.*?)"
        rf"end\s+subroutine\s+{re.escape(name)}",
        flags=re.IGNORECASE | re.DOTALL,
    )
    match = pattern.search(source)
    if not match:
        raise RuntimeError(f"{SOURCE_PATH} is missing subroutine {name}")
    return match.group("body")


def source_keys() -> dict[str, list[str]]:
    source = SOURCE_PATH.read_text(encoding="utf-8")
    groups: dict[str, list[str]] = {}
    for block, subroutine in BLOCK_SUBROUTINES.items():
        keys = sorted({match.group("key").upper() for match in CASE_KEY_RE.finditer(
            subroutine_text(source, subroutine)
        )})
        groups[block] = keys

    flow_cond = set(groups["FLOW_COND"])
    hit_keys = sorted(flow_cond & HIT_KEYS)
    groups["FLOW_COND"] = sorted(flow_cond - HIT_KEYS)
    groups["HIT"] = hit_keys
    return groups


def documented_keys() -> dict[str, list[str]]:
    groups: dict[str, list[str]] = {}
    in_section = False
    for line in DOC_PATH.read_text(encoding="utf-8").splitlines():
        if line.strip() == "## Parser Key Inventory":
            in_section = True
            continue
        if in_section and line.startswith("## "):
            break
        if not in_section:
            continue
        match = ROW_RE.match(line)
        if not match:
            continue
        block = match.group("block")
        keys = sorted({key_match.group("key") for key_match in KEY_RE.finditer(
            match.group("keys")
        )})
        groups[block] = keys
    return groups


def markdown_rows(groups: dict[str, list[str]]) -> str:
    lines = [
        "| Parser block | Keys parsed from `input_util.f90` |",
        "| --- | --- |",
    ]
    for block in [
        "DOMAIN",
        "MODEL",
        "CORIOLIS",
        "FLOW_COND",
        "OUTPUT",
        "TURBINES",
        "SCALARS",
        "HIT",
    ]:
        keys = ", ".join(f"`{key}`" for key in groups[block])
        lines.append(f"| `{block}` | {keys} |")
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Check GPU coverage audit lesgo.conf parser-key coverage."
    )
    parser.add_argument(
        "--list",
        action="store_true",
        help="print the expected Markdown parser-key inventory table",
    )
    args = parser.parse_args()

    expected = source_keys()
    if args.list:
        print(markdown_rows(expected))
        return 0

    documented = documented_keys()
    ok = True

    missing_blocks = sorted(set(expected) - set(documented))
    stale_blocks = sorted(set(documented) - set(expected))
    if missing_blocks:
        print(f"{DOC_PATH} is missing parser-key inventory rows:")
        for block in missing_blocks:
            print(f"  {block}")
        ok = False
    if stale_blocks:
        print(f"{DOC_PATH} has stale parser-key inventory rows:")
        for block in stale_blocks:
            print(f"  {block}")
        ok = False

    for block in sorted(set(expected) & set(documented)):
        expected_keys = set(expected[block])
        documented_keys_set = set(documented[block])
        missing = sorted(expected_keys - documented_keys_set)
        stale = sorted(documented_keys_set - expected_keys)
        if missing:
            print(f"{DOC_PATH} is missing `{block}` keys:")
            for key in missing:
                print(f"  {key}")
            ok = False
        if stale:
            print(f"{DOC_PATH} has stale `{block}` keys:")
            for key in stale:
                print(f"  {key}")
            ok = False

    if not ok:
        print("Run this to regenerate the expected table:")
        print("  python3 tools/check_lesgo_conf_coverage_docs.py --list")
        return 1

    total_keys = sum(len(keys) for keys in expected.values())
    print(
        "lesgo.conf GPU coverage documentation check passed "
        f"({len(expected)} parser groups, {total_keys} keys)."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
