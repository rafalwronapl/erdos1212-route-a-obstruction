#!/usr/bin/env python3
"""Phase-0 audit for Route A gates: omega tails and DP-bank shape.

This intentionally reuses the mortality row-DP transition, but keeps the audit
separate so existing mortality artifacts remain comparable.
"""

from __future__ import annotations

import argparse
import json
import math
from collections import Counter
from datetime import datetime
from pathlib import Path
from typing import Any

from sympy import factorint

from prime_square_segment_a_layered_transition_audit import x_c_for
from prime_square_segment_a_mortality_regeneration_profile import transition_with_deaths
from prime_square_segment_a_phase0_utils import omega_ge5_exact_sieve
from prime_square_segment_b_generic_sentinel_verify import DEFAULT_D
from prime_square_segment_b_row_dp import good


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def contiguous_runs(sorted_values: list[int]) -> list[list[int]]:
    if not sorted_values:
        return []
    runs: list[list[int]] = [[sorted_values[0]]]
    for value in sorted_values[1:]:
        if value == runs[-1][-1] + 1:
            runs[-1].append(value)
        else:
            runs.append([value])
    return runs


def shape_stats(good_bs: list[int], selected_bs: list[int]) -> dict[str, Any]:
    selected = sorted(selected_bs)
    if not selected:
        return {
            "count": 0,
            "maxrun_index": 0,
            "maxrun_B": 0,
            "gapratio_index": 0.0,
            "gapratio_B": 0.0,
            "fragment_count_index": 0,
            "fragment_count_B": 0,
        }

    index_by_b = {b: i for i, b in enumerate(good_bs)}
    selected_indexes = sorted(index_by_b[b] for b in selected)
    index_runs = contiguous_runs(selected_indexes)

    b_step = math.gcd(*(abs(good_bs[i + 1] - good_bs[i]) for i in range(len(good_bs) - 1))) if len(good_bs) > 1 else 1
    normalized_bs = sorted((b - good_bs[0]) // b_step for b in selected)
    b_runs = contiguous_runs(normalized_bs)

    index_span = selected_indexes[-1] - selected_indexes[0] + 1
    b_span = normalized_bs[-1] - normalized_bs[0] + 1
    return {
        "count": len(selected),
        "maxrun_index": max(len(run) for run in index_runs),
        "maxrun_B": max(len(run) for run in b_runs),
        "gapratio_index": len(selected) / index_span,
        "gapratio_B": len(selected) / b_span,
        "fragment_count_index": len(index_runs),
        "fragment_count_B": len(b_runs),
    }


def update_floor_record(current: dict[str, Any] | None, y: int, row_offset: int, stats: dict[str, Any]) -> dict[str, Any]:
    record = {"y": y, "row_offset": row_offset, **stats}
    if current is None:
        return record
    if stats["maxrun_index"] < current["maxrun_index"]:
        return record
    if stats["maxrun_index"] == current["maxrun_index"] and stats["gapratio_index"] < current["gapratio_index"]:
        return record
    return current


def omega_record(y: int) -> dict[str, Any]:
    factors = sorted(int(p) for p in factorint(y + 1))
    ge5 = [p for p in factors if p >= 5]
    return {
        "y": y,
        "n": y + 1,
        "n_factors": factors,
        "omega": len(factors),
        "omega_ge5": len(ge5),
        "n_factors_ge5": ge5,
    }


def audit(
    *,
    p: int,
    q: int,
    b_min: int,
    b_max: int,
    alpha: int,
    radius: int,
    width: int,
    bank_set: set[int],
    dstar: set[int],
    force_composite_xc: bool,
    omega_stride: int,
    omega_exact: bool,
    max_rows: int | None,
    progress_json: Path | None,
) -> dict[str, Any]:
    y_start = p * p
    y_end = q * q
    rows_total = y_end - y_start
    rows_to_scan = rows_total if max_rows is None else min(rows_total, max_rows)
    good_bs = [b for b in range(b_min, b_max + 1) if b % 2 == 0 and b % 3 == p % 3]

    x_meta: dict[int, int] = {}
    states: dict[int, set[int]] = {}
    for b in good_bs:
        _raw, x_c = x_c_for(p, b, alpha, radius, force_composite_xc)
        x_meta[b] = x_c
        start = {d for d in bank_set if good(x_c + d, y_start)}
        if start:
            states[b] = start

    max_omega: dict[str, Any] | None = None
    max_omega_ge5: dict[str, Any] | None = None
    omega_hist: Counter[int] = Counter()
    omega_ge5_hist: Counter[int] = Counter()
    omega_rows_sampled = 0
    exact_omega_ge5: dict[str, Any] | None = None
    if omega_exact:
        exact_omega_ge5 = omega_ge5_exact_sieve(y_start, y_end, q)

    min_live_shape: dict[str, Any] | None = None
    min_dstar_shape: dict[str, Any] | None = None
    max_live_fragments: dict[str, Any] | None = None
    max_dstar_fragments: dict[str, Any] | None = None
    min_living_dstar_pairs: int | None = None
    death_records: dict[int, dict[str, Any]] = {}

    def observe(y: int) -> None:
        nonlocal max_omega, max_omega_ge5, omega_rows_sampled
        nonlocal min_live_shape, min_dstar_shape, max_live_fragments, max_dstar_fragments
        nonlocal min_living_dstar_pairs

        row_offset = y - y_start
        if not omega_exact and row_offset % omega_stride == 0:
            rec = omega_record(y)
            omega_rows_sampled += 1
            omega_hist[rec["omega"]] += 1
            omega_ge5_hist[rec["omega_ge5"]] += 1
            if max_omega is None or rec["omega"] > max_omega["omega"]:
                max_omega = rec
            if max_omega_ge5 is None or rec["omega_ge5"] > max_omega_ge5["omega_ge5"]:
                max_omega_ge5 = rec

        live_bs = [b for b in good_bs if b in states]
        dstar_bs = [b for b in live_bs if states[b] & dstar]
        living_dstar_pairs = sum(len(states[b] & dstar) for b in live_bs)
        min_living_dstar_pairs = (
            living_dstar_pairs
            if min_living_dstar_pairs is None
            else min(min_living_dstar_pairs, living_dstar_pairs)
        )

        live_stats = shape_stats(good_bs, live_bs)
        dstar_stats = shape_stats(good_bs, dstar_bs)
        min_live_shape = update_floor_record(min_live_shape, y, row_offset, live_stats)
        min_dstar_shape = update_floor_record(min_dstar_shape, y, row_offset, dstar_stats)

        live_frag_record = {"y": y, "row_offset": row_offset, **live_stats}
        dstar_frag_record = {"y": y, "row_offset": row_offset, **dstar_stats}
        if max_live_fragments is None or live_stats["fragment_count_index"] > max_live_fragments["fragment_count_index"]:
            max_live_fragments = live_frag_record
        if max_dstar_fragments is None or dstar_stats["fragment_count_index"] > max_dstar_fragments["fragment_count_index"]:
            max_dstar_fragments = dstar_frag_record

    observe(y_start)
    for row_index, y in enumerate(range(y_start, y_start + rows_to_scan), start=1):
        before = states
        states, deaths = transition_with_deaths(states, x_meta, y, width, bank_set)
        for b, rec in deaths.items():
            if b not in death_records:
                rec["Dstar_before_count"] = len(before.get(b, set()) & dstar)
                death_records[b] = rec
        observe(y + 1)
        if not states:
            break
        if progress_json and (row_index == 1 or row_index % 25000 == 0 or row_index == rows_to_scan):
            write_json(
                progress_json,
                {
                    "status": "running",
                    "generated_at": datetime.now().isoformat(timespec="seconds"),
                    "p": p,
                    "q": q,
                    "rows_completed": row_index,
                    "rows_to_scan": rows_to_scan,
                    "omega_rows_sampled": omega_rows_sampled,
                    "death_count": len(death_records),
                    "max_omega": max_omega,
                    "max_omega_ge5": max_omega_ge5,
                    "min_live_shape": min_live_shape,
                    "min_dstar_shape": min_dstar_shape,
                },
            )

    logp2 = math.log(p) ** 2
    kill_omega = math.floor(math.log2(logp2)) + 1
    q1_record = exact_omega_ge5["max_omega_ge5"] if exact_omega_ge5 else max_omega_ge5
    pass_omega_ge5 = q1_record is not None and q1_record["omega_ge5"] < kill_omega

    payload = {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "p": p,
        "q": q,
        "rows_total_full_interval": rows_total,
        "rows_scanned": rows_to_scan,
        "scan_complete_interval": rows_to_scan == rows_total,
        "B_range": [b_min, b_max],
        "good_class_B": good_bs,
        "radius": radius,
        "width": width,
        "bank_set": sorted(bank_set),
        "Dstar": sorted(dstar),
        "omega_stride": omega_stride,
        "omega_exact": omega_exact,
        "omega_rows_sampled": omega_rows_sampled,
        "logp_squared": logp2,
        "kill_line_omega_ge5": kill_omega,
        "phase0_Q1_pass_on_sample": pass_omega_ge5,
        "phase0_Q1_exact": exact_omega_ge5,
        "max_omega": max_omega,
        "max_omega_ge5": q1_record,
        "omega_hist": {str(k): v for k, v in sorted(omega_hist.items())},
        "omega_ge5_hist": {str(k): v for k, v in sorted(omega_ge5_hist.items())},
        "death_count": len(death_records),
        "surviving_B_count": len([b for b in good_bs if b not in death_records]),
        "death_records": {str(k): v for k, v in sorted(death_records.items())},
        "min_living_Dstar_pairs": min_living_dstar_pairs,
        "min_live_shape": min_live_shape,
        "min_dstar_shape": min_dstar_shape,
        "max_live_fragments": max_live_fragments,
        "max_dstar_fragments": max_dstar_fragments,
        "phase0_Q2_live_maxrun_ratio": (
            None if min_live_shape is None or len(good_bs) == 0 else min_live_shape["maxrun_index"] / len(good_bs)
        ),
        "phase0_Q2_dstar_maxrun_ratio": (
            None if min_dstar_shape is None or len(good_bs) == 0 else min_dstar_shape["maxrun_index"] / len(good_bs)
        ),
    }
    return payload


def write_md(path: Path, payload: dict[str, Any]) -> None:
    q1 = payload["max_omega_ge5"] or {}
    lines = [
        "# Segment A Phase-0 Route-A Audit",
        "",
        f"- generated_at: `{payload['generated_at']}`",
        f"- p: `{payload['p']}`",
        f"- q: `{payload['q']}`",
        f"- rows_scanned: `{payload['rows_scanned']}` / `{payload['rows_total_full_interval']}`",
        f"- scan_complete_interval: `{payload['scan_complete_interval']}`",
        f"- omega_stride: `{payload['omega_stride']}`",
        f"- omega_exact: `{payload['omega_exact']}`",
        f"- radius: `{payload['radius']}`",
        f"- good_class_B_count: `{len(payload['good_class_B'])}`",
        "",
        "## Q1 Omega Gate",
        "",
        f"- log(p)^2: `{payload['logp_squared']:.6f}`",
        f"- kill_line_omega_ge5: `{payload['kill_line_omega_ge5']}`",
        f"- max_omega_ge5: `{q1.get('omega_ge5')}` at y `{q1.get('y')}`",
        f"- max_omega_ge5 factors: `{q1.get('n_factors_ge5')}`",
        f"- phase0_Q1_pass_on_sample: `{payload['phase0_Q1_pass_on_sample']}`",
        "",
        "## Q2 Shape Gate",
        "",
        "Live-B floor:",
        "",
        "```json",
        json.dumps(payload["min_live_shape"], indent=2, sort_keys=True),
        "```",
        "",
        "Dstar-nonempty-B floor:",
        "",
        "```json",
        json.dumps(payload["min_dstar_shape"], indent=2, sort_keys=True),
        "```",
        "",
        "## Mortality",
        "",
        f"- death_count: `{payload['death_count']}`",
        f"- surviving_B_count: `{payload['surviving_B_count']}`",
        f"- min_living_Dstar_pairs: `{payload['min_living_Dstar_pairs']}`",
        "",
        "## Histograms",
        "",
        "```json",
        json.dumps(
            {
                "omega_hist": payload["omega_hist"],
                "omega_ge5_hist": payload["omega_ge5_hist"],
            },
            indent=2,
            sort_keys=True,
        ),
        "```",
        "",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p", type=int, required=True)
    parser.add_argument("--q", type=int, required=True)
    parser.add_argument("--b-min", type=int, default=44)
    parser.add_argument("--b-max", type=int, default=104)
    parser.add_argument("--alpha", type=int, default=86)
    parser.add_argument("--radius", type=int, default=50)
    parser.add_argument("--width", type=int, default=30)
    parser.add_argument("--force-composite-xc", action="store_true")
    parser.add_argument("--dstar", default="-26,-22,-20,-10,-8,4,8,10,14")
    parser.add_argument("--bank-even-radius", type=int)
    parser.add_argument("--dstar-even-radius", type=int)
    parser.add_argument("--omega-stride", type=int, default=1)
    parser.add_argument("--omega-exact", action="store_true")
    parser.add_argument("--max-rows", type=int)
    parser.add_argument("--progress-json")
    parser.add_argument("--output-json", required=True)
    parser.add_argument("--output-md", required=True)
    args = parser.parse_args()

    if args.omega_stride <= 0:
        raise SystemExit("--omega-stride must be positive")
    if args.bank_even_radius is not None:
        bank_set = {d for d in range(-args.bank_even_radius, args.bank_even_radius + 1, 2)}
    else:
        bank_set = set(DEFAULT_D)
    if args.dstar_even_radius is not None:
        dstar = {d for d in range(-args.dstar_even_radius, args.dstar_even_radius + 1, 2)}
    else:
        dstar = {int(item) for item in args.dstar.split(",") if item.strip()}
    if not dstar <= bank_set:
        raise SystemExit("Dstar must be a subset of the bank set")
    if max(abs(d) for d in bank_set) > args.width:
        raise SystemExit("--width must cover the bank set")

    payload = audit(
        p=args.p,
        q=args.q,
        b_min=args.b_min,
        b_max=args.b_max,
        alpha=args.alpha,
        radius=args.radius,
        width=args.width,
        bank_set=bank_set,
        dstar=dstar,
        force_composite_xc=args.force_composite_xc,
        omega_stride=args.omega_stride,
        omega_exact=args.omega_exact,
        max_rows=args.max_rows,
        progress_json=Path(args.progress_json) if args.progress_json else None,
    )
    write_json(Path(args.output_json), payload)
    write_md(Path(args.output_md), payload)
    if args.progress_json:
        write_json(Path(args.progress_json), {**payload, "status": "complete", "percent": 100.0})
    print(f"WROTE {args.output_json}")
    print(f"WROTE {args.output_md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
