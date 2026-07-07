#!/usr/bin/env python3
"""Import archived p0 public-case CPU/GPU runs into the evidence ledger."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

from repo_paths import ROOT, repo_path

DEFAULT_EVIDENCE = repo_path("docs", "gpu_validation_evidence.json")
DEFAULT_MANIFEST = repo_path("docs", "gpu_benchmark_manifest.json")


@dataclass(frozen=True)
class P0Case:
    selector: str
    case_id: str
    variant: str | None
    case_dir: Path
    cpu_label: str
    gpu_label: str
    grid: str
    requires_module_check: bool


P0_CASES: dict[str, P0Case] = {
    "les_core_channel": P0Case(
        selector="les_core_channel",
        case_id="les_core_channel",
        variant=None,
        case_dir=Path("test-cases/channel_flow"),
        cpu_label="baseline_cpu",
        gpu_label="baseline_gpu",
        grid="240^3",
        requires_module_check=False,
    ),
    "adm_disk": P0Case(
        selector="adm_disk",
        case_id="adm_disk",
        variant=None,
        case_dir=Path("test-cases/adm_disk"),
        cpu_label="baseline_cpu",
        gpu_label="baseline_gpu",
        grid="240^3",
        requires_module_check=True,
    ),
    "atm_line.structure_off": P0Case(
        selector="atm_line.structure_off",
        case_id="atm_line",
        variant="structure_off",
        case_dir=Path("test-cases/atm_line"),
        cpu_label="structure_off_cpu",
        gpu_label="structure_off_gpu",
        grid="240^3",
        requires_module_check=True,
    ),
    "atm_line.structure_on": P0Case(
        selector="atm_line.structure_on",
        case_id="atm_line",
        variant="structure_on",
        case_dir=Path("test-cases/atm_line"),
        cpu_label="structure_on_cpu",
        gpu_label="structure_on_gpu",
        grid="240^3",
        requires_module_check=True,
    ),
}


def archive_log(root: Path, case: P0Case, label: str) -> Path:
    return root / case.case_dir / "run-archives" / label / f"lesgo_{label}.log"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def manifest_evidence_items(manifest: dict, case_id: str) -> list[str]:
    for row in manifest.get("cases", []):
        if isinstance(row, dict) and row.get("id") == case_id:
            return [
                item
                for item in row.get("required_evidence", [])
                if isinstance(item, str) and item
            ]
    raise KeyError(f"manifest has no case id `{case_id}`")


def parse_module_check(value: str) -> tuple[str, str]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("module checks must use KEY=TEXT")
    key, text = value.split("=", 1)
    key = key.strip()
    text = text.strip()
    if not key or not text:
        raise argparse.ArgumentTypeError("module checks must use non-empty KEY=TEXT")
    return key, text


def build_pair_import_command(
    *,
    case: P0Case,
    evidence: Path,
    manifest: dict,
    archive_root: Path,
    source: str,
    module_checks: dict[str, str],
    public_record: str,
    dry_run: bool,
) -> list[str]:
    cpu_log = archive_log(archive_root, case, case.cpu_label)
    gpu_log = archive_log(archive_root, case, case.gpu_label)
    missing = [str(path) for path in [cpu_log, gpu_log] if not path.exists()]
    if missing:
        raise FileNotFoundError(
            f"{case.selector} is missing archived log(s): " + ", ".join(missing)
        )

    command = [
        sys.executable,
        str(ROOT / "tools/import_lesgo_timing_pair.py"),
        "--evidence",
        str(evidence),
        "--case-id",
        case.case_id,
    ]
    if case.variant:
        command.extend(["--variant", case.variant])
    command.extend(
        [
            "--cpu-log",
            str(cpu_log),
            "--gpu-log",
            str(gpu_log),
            "--grid",
            case.grid,
            "--source",
            source,
            "--compare-diagnostics",
            "--public-record",
            public_record,
        ]
    )

    if case.requires_module_check:
        check = module_checks.get(case.selector) or module_checks.get(case.case_id)
        if not check:
            raise ValueError(
                f"{case.selector} requires --module-check {case.selector}=TEXT "
                "with turbine/module-specific correctness evidence"
            )
        command.extend(["--correctness-check", check])

    for item in manifest_evidence_items(manifest, case.case_id):
        command.extend(["--evidence-item", item])
    if dry_run:
        command.append("--dry-run")
    return command


def shell_quote(value: str) -> str:
    if not value:
        return "''"
    if all(ch.isalnum() or ch in "-_./:=+^" for ch in value):
        return value
    return "'" + value.replace("'", "'\"'\"'") + "'"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Import archived public p0 CPU/GPU runs from run-archives into "
            "docs/gpu_validation_evidence.json."
        )
    )
    parser.add_argument("--evidence", type=Path, default=DEFAULT_EVIDENCE)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument(
        "--archive-root",
        type=Path,
        default=ROOT,
        help="repository root containing test-cases/*/run-archives",
    )
    parser.add_argument("--source", required=True)
    parser.add_argument(
        "--case",
        action="append",
        choices=sorted(P0_CASES),
        help="p0 case selector; repeatable; defaults to all p0 selectors",
    )
    parser.add_argument(
        "--module-check",
        action="append",
        type=parse_module_check,
        default=[],
        metavar="KEY=TEXT",
        help=(
            "module-specific correctness check for ADM/ATM selectors, for "
            "example adm_disk='turbine output matched CPU reference'"
        ),
    )
    parser.add_argument(
        "--public-record",
        choices=["true", "false"],
        default="false",
        help="value passed through to the evidence ledger",
    )
    parser.add_argument(
        "--skip-missing",
        action="store_true",
        help="skip selected cases whose archived logs are not present",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="print import commands and paired-import dry-run output only",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    manifest = load_json(args.manifest)
    archive_root = args.archive_root.resolve()
    selectors = args.case or sorted(P0_CASES)
    module_checks = dict(args.module_check)
    ran = 0
    skipped = 0

    for selector in selectors:
        case = P0_CASES[selector]
        try:
            command = build_pair_import_command(
                case=case,
                evidence=args.evidence,
                manifest=manifest,
                archive_root=archive_root,
                source=args.source,
                module_checks=module_checks,
                public_record=args.public_record,
                dry_run=args.dry_run,
            )
        except FileNotFoundError as exc:
            if args.skip_missing:
                print(f"SKIP: {exc}")
                skipped += 1
                continue
            print(exc)
            return 1
        except (KeyError, ValueError) as exc:
            print(exc)
            return 1

        if args.dry_run:
            print(" ".join(shell_quote(part) for part in command))
        result = subprocess.run(
            command,
            cwd=ROOT,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )
        if result.stdout:
            print(result.stdout.rstrip())
        if result.returncode != 0:
            print(f"FAILED: {selector} import exited with {result.returncode}")
            return result.returncode
        ran += 1

    print(f"Imported p0 archived evidence for {ran} selector(s); skipped {skipped}.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
