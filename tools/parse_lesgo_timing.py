#!/usr/bin/env python3
"""Parse LESGO timing and diagnostic fields from stdout/stderr logs."""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path


FLOAT = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?"
PATTERNS = {
    "iteration": re.compile(r"^\s*Iteration:\s*(?P<value>\d+)\s*$", re.IGNORECASE),
    "dt": re.compile(r"^\s*Time step:\s*(?P<value>{})\s*$".format(FLOAT), re.IGNORECASE),
    "dimensional_time": re.compile(
        r"^\s*Dimensional time:\s*(?P<value>{})\s*$".format(FLOAT),
        re.IGNORECASE,
    ),
    "cfl": re.compile(r"^\s*CFL:\s*(?P<value>{})\s*$".format(FLOAT), re.IGNORECASE),
    "divergence": re.compile(
        r"^\s*Velocity divergence metric:\s*(?P<value>{})\s*$".format(FLOAT),
        re.IGNORECASE,
    ),
    "kinetic_energy": re.compile(
        r"^\s*Kinetic energy:\s*(?P<value>{})\s*$".format(FLOAT),
        re.IGNORECASE,
    ),
    "iteration_wall": re.compile(
        r"^\s*Iteration:\s*(?P<value>{})\s*$".format(FLOAT),
        re.IGNORECASE,
    ),
    "cumulative_wall": re.compile(
        r"^\s*Cumulative:\s*(?P<value>{})\s*$".format(FLOAT),
        re.IGNORECASE,
    ),
    "forcing_wall": re.compile(
        r"^\s*Forcing:\s*(?P<value>{})\s*$".format(FLOAT),
        re.IGNORECASE,
    ),
    "cumulative_forcing": re.compile(
        r"^\s*Cumulative Forcing:\s*(?P<value>{})\s*$".format(FLOAT),
        re.IGNORECASE,
    ),
    "final_wall": re.compile(
        r"^\s*Simulation wall time \(s\)\s*:\s*(?P<value>{})\s*$".format(FLOAT),
        re.IGNORECASE,
    ),
    "final_cpu": re.compile(
        r"^\s*Simulation cpu time \(s\)\s*:\s*(?P<value>{})\s*$".format(FLOAT),
        re.IGNORECASE,
    ),
    "atm_measured_subtotal": re.compile(
        r"^\s*ATM measured subtotal:\s*(?P<value>{})\s*$".format(FLOAT),
        re.IGNORECASE,
    ),
    "mpi_exit_status": re.compile(
        r"\bMPI_EXIT_STATUS\s*[:=]\s*(?P<value>-?\d+)\b",
        re.IGNORECASE,
    ),
}

COMPONENT_RE = re.compile(
    r"^\s*(?P<name>Derivatives|Derivatives xy/filter|Derivatives z|"
    r"SGS & Stresses|SGS model/stress build|SGS tzz halo|"
    r"SGS divstress_uv|SGS divstress_w|Convection|Pressure Solver|"
    r"Projection|Other):\s*(?P<value>{})\s*$".format(FLOAT),
    re.IGNORECASE,
)


def read_text(path: str) -> str:
    if path == "-":
        return sys.stdin.read()
    return Path(path).read_text(encoding="utf-8", errors="replace")


def parse_float(value: str) -> float:
    return float(value.replace("D", "E").replace("d", "e"))


def parse_log(text: str, nsteps: int | None) -> dict:
    result: dict[str, object] = {
        "last_iteration": None,
        "last_iteration_wall_s": None,
        "cumulative_wall_s": None,
        "cumulative_average_s_per_step": None,
        "final_wall_s": None,
        "final_cpu_s": None,
        "mpi_exit_status": None,
        "diagnostics": {},
        "component_cumulative_s": {},
    }
    in_wall_times = False

    for raw_line in text.splitlines():
        line = raw_line.strip()
        if line == "Simulation wall times (s):":
            in_wall_times = True
            continue
        if line.startswith("Sub-component Cumulative Times"):
            in_wall_times = False

        match = PATTERNS["iteration"].match(line)
        if match and not in_wall_times:
            result["last_iteration"] = int(match.group("value"))
            continue

        for key in [
            "dt",
            "dimensional_time",
            "cfl",
            "divergence",
            "kinetic_energy",
        ]:
            match = PATTERNS[key].match(line)
            if match:
                diagnostics = result["diagnostics"]
                assert isinstance(diagnostics, dict)
                diagnostics[key] = parse_float(match.group("value"))
                break
        else:
            match = PATTERNS["iteration_wall"].match(line)
            if match and in_wall_times:
                result["last_iteration_wall_s"] = parse_float(match.group("value"))
                continue
            match = PATTERNS["cumulative_wall"].match(line)
            if match and in_wall_times:
                result["cumulative_wall_s"] = parse_float(match.group("value"))
                continue
            match = PATTERNS["forcing_wall"].match(line)
            if match and in_wall_times:
                result["forcing_wall_s"] = parse_float(match.group("value"))
                continue
            match = PATTERNS["cumulative_forcing"].match(line)
            if match and in_wall_times:
                result["cumulative_forcing_s"] = parse_float(match.group("value"))
                continue
            match = PATTERNS["final_wall"].match(line)
            if match:
                result["final_wall_s"] = parse_float(match.group("value"))
                continue
            match = PATTERNS["final_cpu"].match(line)
            if match:
                result["final_cpu_s"] = parse_float(match.group("value"))
                continue
            match = PATTERNS["atm_measured_subtotal"].match(line)
            if match:
                result["atm_measured_subtotal_s"] = parse_float(match.group("value"))
                continue
            match = PATTERNS["mpi_exit_status"].search(line)
            if match:
                result["mpi_exit_status"] = int(match.group("value"))
                continue
            match = COMPONENT_RE.match(line)
            if match:
                components = result["component_cumulative_s"]
                assert isinstance(components, dict)
                name = match.group("name").lower().replace(" ", "_").replace("/", "_")
                name = name.replace("&", "and")
                components[name] = parse_float(match.group("value"))

    steps = nsteps if nsteps is not None else result["last_iteration"]
    cumulative = result["cumulative_wall_s"]
    if isinstance(steps, int) and steps > 0 and isinstance(cumulative, float):
        result["cumulative_average_s_per_step"] = cumulative / steps
    return result


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Parse LESGO timing and diagnostic fields from a log."
    )
    parser.add_argument("log", help="log file path, or '-' for stdin")
    parser.add_argument(
        "--nsteps",
        type=int,
        help="steps to use when computing cumulative average; defaults to last iteration",
    )
    parser.add_argument("--pretty", action="store_true", help="indent JSON output")
    args = parser.parse_args()

    if args.nsteps is not None and args.nsteps <= 0:
        parser.error("--nsteps must be positive")

    result = parse_log(read_text(args.log), args.nsteps)
    print(json.dumps(result, indent=2 if args.pretty else None, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
