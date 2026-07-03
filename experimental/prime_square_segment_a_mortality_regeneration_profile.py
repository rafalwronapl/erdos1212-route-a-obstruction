#!/usr/bin/env python3
"""Profile Segment A lane mortality separately from alive D* regeneration."""

from __future__ import annotations

import argparse
import json
import math
from collections import Counter
from datetime import datetime
from pathlib import Path
from typing import Any

from sympy import factorint

from prime_square_segment_a_layered_transition_audit import component_id_map, x_c_for
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
            "gapratio_index": 0.0,
            "fragment_count_index": 0,
        }
    index_by_b = {b: i for i, b in enumerate(good_bs)}
    selected_indexes = sorted(index_by_b[b] for b in selected)
    runs = contiguous_runs(selected_indexes)
    span = selected_indexes[-1] - selected_indexes[0] + 1
    return {
        "count": len(selected),
        "maxrun_index": max(len(run) for run in runs),
        "gapratio_index": len(selected) / span,
        "fragment_count_index": len(runs),
    }


def update_shape_floor(current: dict[str, Any] | None, y: int, row_offset: int, stats: dict[str, Any]) -> dict[str, Any]:
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
    factors_ge5 = [p for p in factors if p >= 5]
    return {
        "y": y,
        "n": y + 1,
        "n_factors": factors,
        "omega": len(factors),
        "omega_ge5": len(factors_ge5),
        "n_factors_ge5": factors_ge5,
    }


def transition_with_deaths(
    states: dict[int, set[int]],
    x_meta: dict[int, int],
    y: int,
    width: int,
    bank_set: set[int],
) -> tuple[dict[int, set[int]], dict[int, dict[str, Any]]]:
    next_states: dict[int, set[int]] = {}
    deaths: dict[int, dict[str, Any]] = {}
    for b, bank in states.items():
        x_c = x_meta[b]
        good_y = {d for d in range(-width, width + 1) if good(x_c + d, y)}
        comp = component_id_map(good_y, -width, width)
        left_components = {comp[d] for d in bank if d in comp}
        good_next = {d for d in bank_set if good(x_c + d, y + 1)}
        next_bank = {e for e in good_next if e in comp and comp[e] in left_components}
        if next_bank:
            next_states[b] = next_bank
        else:
            deaths[b] = {
                "death_y": y + 1,
                "transition_y": y,
                "n": y + 1,
                "n_factors": sorted(int(p) for p in factorint(y + 1)),
                "bank_before": sorted(bank),
                "Dstar_before_count": None,
                "good_y_count": len(good_y),
                "good_next_count": len(good_next),
                "left_component_count": len(left_components),
            }
    return next_states, deaths


def finish_alive_zero_window(
    windows: list[dict[str, Any]],
    hist: Counter[int],
    b: int,
    start_y: int,
    end_y: int,
    max_windows: int,
) -> None:
    if start_y <= end_y:
        length = end_y - start_y + 1
        hist[length] += 1
        if len(windows) < max_windows:
            windows.append(
                {
                    "B": b,
                    "start_y": start_y,
                    "end_y": end_y,
                    "length": length,
                    "start_n_factors": sorted(int(p) for p in factorint(start_y + 1)),
                    "end_n_factors": sorted(int(p) for p in factorint(end_y + 1)),
                }
            )


def profile(
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
    max_windows: int,
    probe_offsets: set[int],
    probe_half_window: int,
    progress_path: Path | None,
    phase0_audit: bool,
    omega_stride: int,
    omega_exact: bool,
) -> dict[str, Any]:
    y_start = p * p
    y_end = q * q
    rows_total = y_end - y_start
    good_bs = [b for b in range(b_min, b_max + 1) if b % 2 == 0 and b % 3 == p % 3]

    x_meta: dict[int, int] = {}
    states: dict[int, set[int]] = {}
    for b in good_bs:
        _raw, x_c = x_c_for(p, b, alpha, radius, force_composite_xc)
        x_meta[b] = x_c
        start = {d for d in DEFAULT_D if good(x_c + d, y_start)}
        if bank_set != set(DEFAULT_D):
            start = {d for d in bank_set if good(x_c + d, y_start)}
        if start:
            states[b] = start

    alive = {b for b in good_bs if b in states}
    death_records: dict[int, dict[str, Any]] = {}
    alive_zero_start: dict[int, int | None] = {b: None for b in good_bs}
    alive_zero_windows: list[dict[str, Any]] = []
    alive_zero_len_hist: Counter[int] = Counter()
    simultaneous_alive_zero_hist: Counter[int] = Counter()
    dead_count_hist: Counter[int] = Counter()
    living_nonempty_b_hist: Counter[int] = Counter()
    max_alive_zero_len_by_b = {b: 0 for b in good_bs}
    total_alive_zero_rows_by_b: Counter[int] = Counter()

    min_living_dstar_pairs: int | None = None
    min_living_nonempty_b: int | None = None
    min_living_record: dict[str, Any] | None = None
    max_simultaneous_alive_zero = 0
    max_omega: dict[str, Any] | None = None
    max_omega_ge5: dict[str, Any] | None = None
    omega_rows_sampled = 0
    omega_hist: Counter[int] = Counter()
    omega_ge5_hist: Counter[int] = Counter()
    exact_omega_ge5: dict[str, Any] | None = None
    if phase0_audit and omega_exact:
        exact_omega_ge5 = omega_ge5_exact_sieve(y_start, y_end, q)
    min_live_shape: dict[str, Any] | None = None
    min_dstar_shape: dict[str, Any] | None = None
    max_live_fragments: dict[str, Any] | None = None
    max_dstar_fragments: dict[str, Any] | None = None
    probe_windows: dict[int, dict[str, Any]] = {
        off: {
            "target_offset": off,
            "target_y": y_start + off,
            "half_window": probe_half_window,
            "center": None,
            "min_living_Dstar_pairs": None,
            "min_living_Dstar_nonempty_B": None,
            "min_Dstar_counts_by_B": {},
            "min_bank_sizes_by_B": {},
        }
        for off in sorted(probe_offsets)
    }

    def observe(y: int, current_states: dict[int, set[int]]) -> None:
        nonlocal min_living_dstar_pairs, min_living_nonempty_b, min_living_record, max_simultaneous_alive_zero
        nonlocal max_omega, max_omega_ge5, omega_rows_sampled
        nonlocal min_live_shape, min_dstar_shape, max_live_fragments, max_dstar_fragments
        row_offset = y - y_start
        live_bs = [b for b in good_bs if b in current_states]
        dead_bs = [b for b in good_bs if b not in current_states]
        counts = {b: len(current_states.get(b, set()) & dstar) for b in live_bs}
        bank_sizes = {b: len(current_states.get(b, set())) for b in live_bs}
        alive_zero_bs = [b for b, count in counts.items() if count == 0]
        living_nonempty_b = len(live_bs) - len(alive_zero_bs)
        living_pairs = sum(counts.values())

        if phase0_audit:
            if not omega_exact and row_offset % omega_stride == 0:
                rec = omega_record(y)
                omega_rows_sampled += 1
                omega_hist[rec["omega"]] += 1
                omega_ge5_hist[rec["omega_ge5"]] += 1
                if max_omega is None or rec["omega"] > max_omega["omega"]:
                    max_omega = rec
                if max_omega_ge5 is None or rec["omega_ge5"] > max_omega_ge5["omega_ge5"]:
                    max_omega_ge5 = rec

            dstar_nonempty_bs = [b for b in live_bs if counts[b] > 0]
            live_stats = shape_stats(good_bs, live_bs)
            dstar_stats = shape_stats(good_bs, dstar_nonempty_bs)
            min_live_shape = update_shape_floor(min_live_shape, y, row_offset, live_stats)
            min_dstar_shape = update_shape_floor(min_dstar_shape, y, row_offset, dstar_stats)
            live_fragment_record = {"y": y, "row_offset": row_offset, **live_stats}
            dstar_fragment_record = {"y": y, "row_offset": row_offset, **dstar_stats}
            if (
                max_live_fragments is None
                or live_stats["fragment_count_index"] > max_live_fragments["fragment_count_index"]
            ):
                max_live_fragments = live_fragment_record
            if (
                max_dstar_fragments is None
                or dstar_stats["fragment_count_index"] > max_dstar_fragments["fragment_count_index"]
            ):
                max_dstar_fragments = dstar_fragment_record

        simultaneous_alive_zero_hist[len(alive_zero_bs)] += 1
        dead_count_hist[len(dead_bs)] += 1
        living_nonempty_b_hist[living_nonempty_b] += 1
        max_simultaneous_alive_zero = max(max_simultaneous_alive_zero, len(alive_zero_bs))

        if (
            min_living_dstar_pairs is None
            or living_pairs < min_living_dstar_pairs
            or living_nonempty_b < (min_living_nonempty_b if min_living_nonempty_b is not None else 10**9)
        ):
            min_living_dstar_pairs = living_pairs if min_living_dstar_pairs is None else min(min_living_dstar_pairs, living_pairs)
            min_living_nonempty_b = (
                living_nonempty_b
                if min_living_nonempty_b is None
                else min(min_living_nonempty_b, living_nonempty_b)
            )
            min_living_record = {
                "y": y,
                "n": y + 1,
                "n_factors": sorted(int(p0) for p0 in factorint(y + 1)),
                "live_B_count": len(live_bs),
                "dead_B_count": len(dead_bs),
                "living_Dstar_nonempty_B": living_nonempty_b,
                "living_Dstar_pairs": living_pairs,
                "alive_zero_B": alive_zero_bs,
                "dead_B": dead_bs,
                "Dstar_counts": {str(b): counts[b] for b in sorted(counts)},
            }

        for b in live_bs:
            if counts[b] == 0:
                total_alive_zero_rows_by_b[b] += 1
                if alive_zero_start[b] is None:
                    alive_zero_start[b] = y
            elif alive_zero_start[b] is not None:
                start = int(alive_zero_start[b])
                length = y - start
                max_alive_zero_len_by_b[b] = max(max_alive_zero_len_by_b[b], length)
                finish_alive_zero_window(
                    alive_zero_windows,
                    alive_zero_len_hist,
                    b,
                    start,
                    y - 1,
                    max_windows,
                )
                alive_zero_start[b] = None

        for target_offset, probe in probe_windows.items():
            if abs(row_offset - target_offset) > probe_half_window:
                continue
            record = {
                "y": y,
                "row_offset": row_offset,
                "n": y + 1,
                "n_factors": sorted(int(p0) for p0 in factorint(y + 1)),
                "living_Dstar_pairs": living_pairs,
                "living_Dstar_nonempty_B": living_nonempty_b,
                "Dstar_counts": {str(b): counts[b] for b in sorted(counts)},
                "bank_sizes": {str(b): bank_sizes[b] for b in sorted(bank_sizes)},
            }
            if row_offset == target_offset:
                probe["center"] = record
            if (
                probe["min_living_Dstar_pairs"] is None
                or living_pairs < probe["min_living_Dstar_pairs"]["living_Dstar_pairs"]
            ):
                probe["min_living_Dstar_pairs"] = record
            if (
                probe["min_living_Dstar_nonempty_B"] is None
                or living_nonempty_b < probe["min_living_Dstar_nonempty_B"]["living_Dstar_nonempty_B"]
            ):
                probe["min_living_Dstar_nonempty_B"] = record
            for b in good_bs:
                b_key = str(b)
                c = counts.get(b, 0)
                if (
                    b_key not in probe["min_Dstar_counts_by_B"]
                    or c < probe["min_Dstar_counts_by_B"][b_key]["Dstar_count"]
                ):
                    probe["min_Dstar_counts_by_B"][b_key] = {
                        "B": b,
                        "Dstar_count": c,
                        "bank_size": bank_sizes.get(b, 0),
                        "y": y,
                        "row_offset": row_offset,
                    }
                s = bank_sizes.get(b, 0)
                if (
                    b_key not in probe["min_bank_sizes_by_B"]
                    or s < probe["min_bank_sizes_by_B"][b_key]["bank_size"]
                ):
                    probe["min_bank_sizes_by_B"][b_key] = {
                        "B": b,
                        "bank_size": s,
                        "Dstar_count": c,
                        "y": y,
                        "row_offset": row_offset,
                    }

    observe(y_start, states)
    for row_index, y in enumerate(range(y_start, y_end), start=1):
        before_states = states
        states, deaths = transition_with_deaths(states, x_meta, y, width, bank_set)
        for b, rec in deaths.items():
            if b not in death_records:
                rec["Dstar_before_count"] = len(before_states.get(b, set()) & dstar)
                death_records[b] = rec
            if alive_zero_start.get(b) is not None:
                start = int(alive_zero_start[b])
                length = y - start
                max_alive_zero_len_by_b[b] = max(max_alive_zero_len_by_b[b], length)
                finish_alive_zero_window(
                    alive_zero_windows,
                    alive_zero_len_hist,
                    b,
                    start,
                    y - 1,
                    max_windows,
                )
                alive_zero_start[b] = None

        observe(y + 1, states)

        if progress_path and (row_index == 1 or row_index % 25000 == 0 or row_index == rows_total):
            progress_payload = {
                "status": "running",
                "generated_at": datetime.now().isoformat(timespec="seconds"),
                "p": p,
                "q": q,
                "rows_completed": row_index,
                "rows_total": rows_total,
                "percent": round(100.0 * row_index / rows_total, 2),
                "dead_B_count": len(death_records),
                "min_living_Dstar_pairs": min_living_dstar_pairs,
                "min_living_Dstar_nonempty_B": min_living_nonempty_b,
                "max_simultaneous_alive_zero_B": max_simultaneous_alive_zero,
            }
            if phase0_audit:
                progress_payload["phase0_audit"] = {
                    "omega_stride": omega_stride,
                    "omega_exact": omega_exact,
                    "omega_rows_sampled": omega_rows_sampled,
                    "max_omega": max_omega,
                    "max_omega_ge5": exact_omega_ge5["max_omega_ge5"] if exact_omega_ge5 else max_omega_ge5,
                    "min_live_shape": min_live_shape,
                    "min_dstar_shape": min_dstar_shape,
                }
            write_json(
                progress_path,
                progress_payload,
            )
        if not states:
            break

    final_y = y_end
    for b, start in alive_zero_start.items():
        if start is not None:
            length = final_y - start + 1
            max_alive_zero_len_by_b[b] = max(max_alive_zero_len_by_b[b], length)
            finish_alive_zero_window(alive_zero_windows, alive_zero_len_hist, b, start, final_y, max_windows)

    payload = {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "p": p,
        "q": q,
        "rows_total": rows_total,
        "B_range": [b_min, b_max],
        "good_class_B": good_bs,
        "Dstar": sorted(dstar),
        "bank_set": sorted(bank_set),
        "death_count": len(death_records),
        "surviving_B_count": len([b for b in good_bs if b not in death_records]),
        "death_records": {str(k): v for k, v in sorted(death_records.items())},
        "min_living_Dstar_pairs": min_living_dstar_pairs,
        "min_living_Dstar_nonempty_B": min_living_nonempty_b,
        "min_living_record": min_living_record,
        "max_simultaneous_alive_zero_B": max_simultaneous_alive_zero,
        "alive_zero_len_hist": {str(k): v for k, v in sorted(alive_zero_len_hist.items())},
        "simultaneous_alive_zero_B_hist": {str(k): v for k, v in sorted(simultaneous_alive_zero_hist.items())},
        "dead_count_hist": {str(k): v for k, v in sorted(dead_count_hist.items())},
        "living_nonempty_B_hist": {str(k): v for k, v in sorted(living_nonempty_b_hist.items())},
        "max_alive_zero_len_by_B": {str(b): max_alive_zero_len_by_b[b] for b in good_bs},
        "total_alive_zero_rows_by_B": {str(b): total_alive_zero_rows_by_b[b] for b in good_bs},
        "sample_alive_zero_windows": alive_zero_windows,
        "probe_windows": probe_windows,
    }
    if phase0_audit:
        logp2 = math.log(p) ** 2
        kill_omega = math.floor(math.log2(logp2)) + 1
        q1_record = exact_omega_ge5["max_omega_ge5"] if exact_omega_ge5 else max_omega_ge5
        payload["phase0_audit"] = {
            "omega_stride": omega_stride,
            "omega_exact": omega_exact,
            "omega_rows_sampled": omega_rows_sampled,
            "logp_squared": logp2,
            "kill_line_omega_ge5": kill_omega,
            "phase0_Q1_pass": q1_record is not None and q1_record["omega_ge5"] < kill_omega,
            "phase0_Q1_exact": exact_omega_ge5,
            "max_omega": max_omega,
            "max_omega_ge5": q1_record,
            "omega_hist": {str(k): v for k, v in sorted(omega_hist.items())},
            "omega_ge5_hist": {str(k): v for k, v in sorted(omega_ge5_hist.items())},
            "min_live_shape": min_live_shape,
            "min_dstar_shape": min_dstar_shape,
            "max_live_fragments": max_live_fragments,
            "max_dstar_fragments": max_dstar_fragments,
            "live_maxrun_ratio": None if min_live_shape is None or not good_bs else min_live_shape["maxrun_index"] / len(good_bs),
            "dstar_maxrun_ratio": None if min_dstar_shape is None or not good_bs else min_dstar_shape["maxrun_index"] / len(good_bs),
        }
    return payload


def write_md(path: Path, payload: dict[str, Any]) -> None:
    lines = [
        "# Segment A Mortality / Regeneration Profile",
        "",
        f"- p: `{payload['p']}`",
        f"- q: `{payload['q']}`",
        f"- rows_total: `{payload['rows_total']}`",
        f"- B_range: `{payload['B_range']}`",
        f"- good_class_B: `{payload['good_class_B']}`",
        "",
        "## Mortality",
        "",
        f"- death_count: `{payload['death_count']}`",
        f"- surviving_B_count: `{payload['surviving_B_count']}`",
        "",
        "## Living-Lane Floors",
        "",
        f"- min_living_Dstar_pairs: `{payload['min_living_Dstar_pairs']}`",
        f"- min_living_Dstar_nonempty_B: `{payload['min_living_Dstar_nonempty_B']}`",
        f"- max_simultaneous_alive_zero_B: `{payload['max_simultaneous_alive_zero_B']}`",
        "",
        "## Death Records",
        "",
        "| B | death_y | n factors | D* before | bank_before |",
        "|---:|---:|---|---:|---|",
    ]
    for b, rec in payload["death_records"].items():
        lines.append(
            f"| {b} | {rec['death_y']} | `{rec['n_factors']}` | {rec['Dstar_before_count']} | `{rec['bank_before']}` |"
        )
    lines.extend(
        [
            "",
            "## Alive D*-Empty Window Lengths",
            "",
            "| length | count |",
            "|---:|---:|",
        ]
    )
    for k, v in payload["alive_zero_len_hist"].items():
        lines.append(f"| {k} | {v} |")
    lines.extend(["", "## Min Living Record", "", "```json", json.dumps(payload["min_living_record"], indent=2), "```", ""])
    if payload.get("probe_windows"):
        lines.extend(["", "## Probe Windows", ""])
        for off, probe in payload["probe_windows"].items():
            lines.extend(
                [
                    f"### offset {off}",
                    "",
                    "```json",
                    json.dumps(probe, indent=2, sort_keys=True),
                    "```",
                    "",
                ]
            )
    if payload.get("phase0_audit"):
        audit = payload["phase0_audit"]
        max_ge5 = audit.get("max_omega_ge5") or {}
        lines.extend(
            [
                "",
                "## Phase-0 Route-A Audit",
                "",
                f"- omega_stride: `{audit['omega_stride']}`",
                f"- omega_exact: `{audit['omega_exact']}`",
                f"- omega_rows_sampled: `{audit['omega_rows_sampled']}`",
                f"- log(p)^2: `{audit['logp_squared']:.6f}`",
                f"- kill_line_omega_ge5: `{audit['kill_line_omega_ge5']}`",
                f"- max_omega_ge5: `{max_ge5.get('omega_ge5')}` at y `{max_ge5.get('y')}`",
                f"- phase0_Q1_pass: `{audit['phase0_Q1_pass']}`",
                "",
                "### Q2 live-bank shape floor",
                "",
                "```json",
                json.dumps(audit["min_live_shape"], indent=2, sort_keys=True),
                "```",
                "",
                "### Q2 Dstar-nonempty shape floor",
                "",
                "```json",
                json.dumps(audit["min_dstar_shape"], indent=2, sort_keys=True),
                "```",
                "",
            ]
        )
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
    parser.add_argument("--max-windows", type=int, default=100)
    parser.add_argument("--probe-offsets", default="")
    parser.add_argument("--probe-half-window", type=int, default=0)
    parser.add_argument("--progress-json")
    parser.add_argument("--phase0-audit", action="store_true")
    parser.add_argument("--omega-stride", type=int, default=1)
    parser.add_argument("--omega-exact", action="store_true")
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
    payload = profile(
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
        max_windows=args.max_windows,
        probe_offsets={int(x) for x in args.probe_offsets.split(",") if x.strip()},
        probe_half_window=args.probe_half_window,
        progress_path=Path(args.progress_json) if args.progress_json else None,
        phase0_audit=args.phase0_audit,
        omega_stride=args.omega_stride,
        omega_exact=args.omega_exact,
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
