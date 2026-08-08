#!/usr/bin/env python3
"""Verify the checked-in Level Set validation matrix cannot silently shrink."""

from __future__ import annotations

from collections import Counter
from pathlib import Path

from prepare_level_set_validation import load_variants


VARIANT_PATH = Path("test-cases/level_set_cubes/validation_variants.json")
REQUIRED_VARIANTS = {
    "sgs_off",
    *(f"sgs_model_{model}" for model in range(1, 6)),
    "model4_beta_off_lbc0",
    "model4_beta_on_lbc1",
    "extrap_tau_log",
    "legacy_extrapolation",
    "modify_dutdn",
    "desired_log_velocity",
    "smooth_3d",
    "sphere_rank_crossing",
    "sphere_rank_crossing_model5",
    "sphere_rank_crossing_4rank",
    "reject_invalid_smooth_mode",
    "reject_smooth_3d_mpi",
    "reject_nonpositive_log_roughness",
    "reject_nonpositive_direct_log_roughness",
    "reject_nonpositive_global_ca_skip",
    "reject_multirank_log_extrapolation",
    "reject_multirank_legacy_extrapolation",
    "reject_multirank_modify_dutdn",
}
STANDARD_PROFILES = {"cpu_mpi", "bridge_mpi", "gpu_mpi_staged"}


def main() -> int:
    variants = load_variants(VARIANT_PATH)
    by_name = {str(item["name"]): item for item in variants}
    if set(by_name) != REQUIRED_VARIANTS:
        missing = sorted(REQUIRED_VARIANTS - set(by_name))
        extra = sorted(set(by_name) - REQUIRED_VARIANTS)
        raise SystemExit(f"Level Set variant set changed: missing={missing}, extra={extra}")

    task_counts = Counter()
    for item in variants:
        task_counts[int(item["nproc"])] += len(item["profiles"])
        settings = item["settings"]
        if settings["dyn_init"] != "1" or settings["cs_count"] != "2":
            raise SystemExit(f"{item['name']}: dynamic SGS update settings changed")

    if dict(task_counts) != {1: 45, 2: 11, 4: 2}:
        raise SystemExit(f"Level Set task counts changed: {dict(task_counts)}")

    for name in ["sgs_off", *(f"sgs_model_{model}" for model in range(1, 6))]:
        if not STANDARD_PROFILES.issubset(set(by_name[name]["profiles"])):
            raise SystemExit(f"{name}: CPU/bridge/full-GPU coverage is incomplete")
    if {"cpu_nompi", "gpu_nompi"} - set(by_name["sgs_model_5"]["profiles"]):
        raise SystemExit("sgs_model_5: non-MPI CPU/GPU coverage is incomplete")
    if {"cpu_nompi", "gpu_nompi"} != set(by_name["smooth_3d"]["profiles"]):
        raise SystemExit("smooth_3d: expected single-rank non-MPI profiles")

    for name in ("sphere_rank_crossing", "sphere_rank_crossing_model5"):
        item = by_name[name]
        if item["geometry"] != "sphere" or int(item["nproc"]) != 2:
            raise SystemExit(f"{name}: rank-crossing geometry contract changed")
        if "gpu_mpi_aware" not in item["profiles"]:
            raise SystemExit(f"{name}: GPU-aware MPI coverage is missing")
    four_rank = by_name["sphere_rank_crossing_4rank"]
    if int(four_rank["nproc"]) != 4 or set(four_rank["profiles"]) != {
        "gpu_mpi_staged", "gpu_mpi_aware"
    }:
        raise SystemExit("sphere_rank_crossing_4rank: staged/aware coverage changed")

    for name, item in by_name.items():
        if name.startswith("reject_") and not item["expected_error"]:
            raise SystemExit(f"{name}: expected rejection text is missing")

    print(
        "Level Set validation-variant check passed "
        f"({len(variants)} variants, {sum(task_counts.values())} tasks)."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
