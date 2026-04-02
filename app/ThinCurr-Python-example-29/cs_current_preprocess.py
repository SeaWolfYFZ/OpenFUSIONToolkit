#!/usr/bin/env python3
"""
Utility for Central Solenoid current waveform.

Features:
1. Check whether the time grid in Central_Solenoid_Current.txt is
   uniformly spaced with 1 ms intervals.
2. Optionally export a compact binary file that is easy to read from
   both C/C++ and Python when coupling to ThinCurr.

Default input layout (tab or space separated, UTF-8 header is ignored):
    Time（ms）Current（A）
    -500    0
    -499    280
    ...

Binary export format (little-endian):
    int64  n_points
    float64 t0_ms      (time of first sample, in ms)
    float64 dt_ms      (uniform time step, in ms)
    float64 I[n_points] (currents in A)
"""

from __future__ import annotations

import argparse
import struct
from pathlib import Path
from typing import List, Tuple


def load_time_current_txt(path: Path) -> Tuple[List[float], List[float]]:
    """Load time (ms) and current (A) from a plain text file."""
    times_ms: List[float] = []
    currents_a: List[float] = []

    with path.open("r", encoding="utf-8") as f:
        for line in f:
            stripped = line.strip()
            if not stripped:
                continue
            # Skip header or comment lines
            if stripped.startswith("#") or "Time" in stripped or "Current" in stripped:
                continue

            parts = stripped.split()
            if len(parts) < 2:
                continue

            try:
                t_ms = float(parts[0])
                I_a = float(parts[1])
            except ValueError:
                # Non-numeric line, skip
                continue

            times_ms.append(t_ms)
            currents_a.append(I_a)

    if not times_ms:
        raise ValueError(f"No valid data lines found in {path}")

    return times_ms, currents_a


def check_uniform_1ms(times_ms: List[float], tol_rel: float = 1e-9) -> Tuple[bool, float, float]:
    """
    Check whether time points are on a uniform grid with 1 ms spacing.

    Returns (is_uniform_1ms, dt_min, dt_max), where dt_min/dt_max are the
    minimal and maximal consecutive differences (in ms).
    """
    if len(times_ms) < 2:
        return True, 0.0, 0.0

    dts = [times_ms[i + 1] - times_ms[i] for i in range(len(times_ms) - 1)]
    dt_min = min(dts)
    dt_max = max(dts)

    # Check uniform grid: all steps essentially equal
    dt_ref = 0.5 * (dt_min + dt_max)
    if dt_ref == 0.0:
        is_uniform = all(abs(dt) < tol_rel for dt in dts)
    else:
        is_uniform = all(abs(dt - dt_ref) <= tol_rel * max(1.0, abs(dt_ref)) for dt in dts)

    # Then check specifically for 1 ms spacing
    is_1ms = abs(dt_ref - 1.0) <= tol_rel * max(1.0, abs(dt_ref))

    return bool(is_uniform and is_1ms), dt_min, dt_max


def export_binary(
    out_path: Path,
    times_ms: List[float],
    currents_a: List[float],
) -> None:
    """
    Export time/current series to a small binary file that is easy to read
    from both C/C++ and Python.

    Layout (little-endian):
        int64   n_points
        float64 t0_ms
        float64 dt_ms
        float64 I[n_points]
    """
    if len(times_ms) != len(currents_a):
        raise ValueError("times and currents must have the same length")
    n = len(times_ms)
    if n == 0:
        raise ValueError("no data to export")
    if n == 1:
        raise ValueError("need at least two samples to define dt_ms")

    dts = [times_ms[i + 1] - times_ms[i] for i in range(n - 1)]
    dt_ms = sum(dts) / float(len(dts))

    header_fmt = "<qdd"  # int64, double, double
    header = struct.pack(header_fmt, n, float(times_ms[0]), float(dt_ms))

    data_fmt = "<" + "d" * n
    data = struct.pack(data_fmt, *[float(I) for I in currents_a])

    with out_path.open("wb") as f:
        f.write(header)
        f.write(data)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Check CS current time grid (1 ms) and optionally export a compact binary waveform file."
    )
    parser.add_argument(
        "--input",
        "-i",
        type=Path,
        default=Path("Central_Solenoid_Current.txt"),
        help="Input time-current text file (default: Central_Solenoid_Current.txt)",
    )
    parser.add_argument(
        "--export-binary",
        "-b",
        type=Path,
        default=None,
        metavar="OUT_BIN",
        help=(
            "If set, export waveform to a simple binary file for C/Python. "
            "Recommended name: Central_Solenoid_Current.bin"
        ),
    )

    args = parser.parse_args()

    times_ms, currents_a = load_time_current_txt(args.input)

    is_uniform_1ms, dt_min, dt_max = check_uniform_1ms(times_ms)

    print(f"Loaded {len(times_ms)} samples from {args.input}")
    print(f"time range (ms): [{min(times_ms)}, {max(times_ms)}]")
    print(f"dt_min (ms): {dt_min}")
    print(f"dt_max (ms): {dt_max}")

    if is_uniform_1ms:
        print("Result: time grid is UNIFORM with Δt ≈ 1 ms.")
    else:
        print("Result: time grid is NOT strictly uniform 1 ms spacing.")

    if args.export_binary is not None:
        out_path: Path = args.export_binary
        export_binary(out_path, times_ms, currents_a)
        print(f"Binary waveform exported to: {out_path}")


if __name__ == "__main__":
    main()

