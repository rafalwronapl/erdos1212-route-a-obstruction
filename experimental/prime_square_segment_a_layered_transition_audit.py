#!/usr/bin/env python3
"""Audit one-step layered transition graphs for a Segment A manifest.

This is a deterministic diagnostic for the Hall/Menger proof route.  For each
row transition y -> y+1 and each good-class B, it builds the bipartite graph
from current reachable offsets d to next offsets e in DEFAULT_D.  An edge exists
when d and e lie in the same good horizontal component in row y and e is good on
row y+1.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from datetime import datetime
from pathlib import Path
from typing import Any

from sympy import isprime

from prime_square_segment_a_family_reservoir_audit import (
    blocker_keys,
    collision_lower_bound,
    exact_independence_number,
)
from prime_square_segment_b_generic_sentinel_verify import DEFAULT_D
from prime_square_segment_b_row_dp import components, good


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def x_c_for(p: int, b: int, alpha: int, radius: int, force_composite: bool) -> tuple[int, int]:
    raw = alpha * p * p + b * p + radius + 1
    x_c = raw
    if force_composite:
        while isprime(x_c):
            x_c += 1
    return raw, x_c


def component_id_map(good_offsets: set[int], min_d: int, max_d: int) -> dict[int, int]:
    mapping: dict[int, int] = {}
    for idx, (left, right) in enumerate(components(good_offsets, min_d, max_d)):
        for d in range(left, right + 1):
            mapping[d] = idx
    return mapping


def audit_manifest(
    manifest: dict[str, Any],
    b_min: int,
    b_max: int,
    alpha: int,
    radius: int,
    width: int,
    force_composite: bool,
    dstar: set[int],
    max_records: int,
    progress_path: Path | None,
    skip_exact_hall: bool,
) -> dict[str, Any]:
    p = int(manifest["p"])
    q = int(manifest["q"])
    y_start = p * p
    y_end = q * q
    bank_set = set(DEFAULT_D)
    good_bs = [b for b in range(b_min, b_max + 1) if b % 2 == 0 and b % 3 == p % 3]

    meta: dict[int, dict[str, int]] = {}
    states: dict[int, set[int]] = {}
    for b in good_bs:
        raw, x_c = x_c_for(p, b, alpha, radius, force_composite)
        meta[b] = {"raw_x_c": raw, "x_c": x_c, "force_shift": x_c - raw}
        start = {d for d in bank_set if good(x_c + d, y_start)}
        if start:
            states[b] = start

    min_right_pairs = sum(len(v) for v in states.values())
    min_dstar_right_pairs = sum(len(v & dstar) for v in states.values())
    min_active_b = len(states)
    min_dstar_active_b = sum(1 for bank in states.values() if bank & dstar)
    min_edge_count = None
    min_record: dict[str, Any] | None = None
    dstar_min_record: dict[str, Any] | None = None
    active_b_hist: Counter[int] = Counter()
    right_pair_hist: Counter[int] = Counter()
    dstar_pair_hist: Counter[int] = Counter()
    thin_records: list[dict[str, Any]] = []
    min_dstar_phi_exact_color_m = None
    min_dstar_phi_record: dict[str, Any] | None = None

    total_rows = y_end - y_start
    for row_index, y in enumerate(range(y_start, y_end), start=1):
        next_states: dict[int, set[int]] = {}
        total_edges = 0
        row_details: list[dict[str, Any]] = []
        dstar_lift_vertices: list[tuple[int, int, int, int]] = []

        for b, bank in states.items():
            x_c = meta[b]["x_c"]
            good_y = {d for d in range(-width, width + 1) if good(x_c + d, y)}
            comp = component_id_map(good_y, -width, width)
            good_next = {d for d in bank_set if good(x_c + d, y + 1)}
            next_bank: set[int] = set()
            edge_count = 0
            left_by_component: Counter[int] = Counter()
            right_by_component: Counter[int] = Counter()

            for d in bank:
                if d in comp:
                    left_by_component[comp[d]] += 1
            for e in good_next:
                if e in comp and left_by_component[comp[e]] > 0:
                    next_bank.add(e)
                    right_by_component[comp[e]] += 1
                    edge_count += left_by_component[comp[e]]
            for e in bank_set & dstar:
                if e in comp and left_by_component[comp[e]] > 0:
                    dstar_lift_vertices.append((len(dstar_lift_vertices), b, e, x_c + e))

            total_edges += edge_count
            if next_bank:
                next_states[b] = next_bank
            row_details.append(
                {
                    "B": b,
                    "left_bank": sorted(bank),
                    "left_size": len(bank),
                    "right_bank": sorted(next_bank),
                    "right_size": len(next_bank),
                    "dstar_right": sorted(next_bank & dstar),
                    "dstar_right_size": len(next_bank & dstar),
                    "transition_edge_count": edge_count,
                    "active_components": len(left_by_component),
                    "right_components": len(right_by_component),
                }
            )

        right_pairs = sum(len(v) for v in next_states.values())
        dstar_right_pairs = sum(len(v & dstar) for v in next_states.values())
        active_b = len(next_states)
        dstar_active_b = sum(1 for bank in next_states.values() if bank & dstar)
        dstar_d_support = {d for bank in next_states.values() for d in bank & dstar}

        active_b_hist[active_b] += 1
        right_pair_hist[right_pairs] += 1
        dstar_pair_hist[dstar_right_pairs] += 1

        row_summary = {
            "y": y,
            "next_y": y + 1,
            "active_B_after": active_b,
            "right_pairs_after": right_pairs,
            "Dstar_active_B_after": dstar_active_b,
            "Dstar_pairs_after": dstar_right_pairs,
            "Dstar_d_support_after": sorted(dstar_d_support),
            "Dstar_d_support_size_after": len(dstar_d_support),
            "Dstar_lift_vertex_count": len(dstar_lift_vertices),
            "transition_edge_count": total_edges,
            "row_details": row_details,
        }

        should_capture_thin_row = (
            active_b <= 4
            or right_pairs <= 15
            or dstar_active_b <= 4
            or dstar_right_pairs <= 13
            or len(dstar_lift_vertices) <= 15
        )
        should_exact_hall = should_capture_thin_row and not skip_exact_hall
        if should_exact_hall:
            collision = collision_lower_bound(dstar_lift_vertices, y + 1)
            alpha_exact = exact_independence_number(dstar_lift_vertices, y + 1)
            if alpha_exact is not None:
                colors: set[str] = set()
                for _idx, _b, _d, x in dstar_lift_vertices:
                    colors.update(blocker_keys(x, y + 1))
                m_color = len(colors)
                collision["alpha_exact"] = alpha_exact
                collision["m_color_exact"] = m_color
                collision["phi_alpha_exact_script_m"] = alpha_exact - collision["m_n"]
                collision["phi_alpha_exact_color_m"] = alpha_exact - m_color
                if (
                    min_dstar_phi_exact_color_m is None
                    or collision["phi_alpha_exact_color_m"] < min_dstar_phi_exact_color_m
                ):
                    min_dstar_phi_exact_color_m = collision["phi_alpha_exact_color_m"]
                    min_dstar_phi_record = {**row_summary, "Dstar_collision": collision}
            row_summary["Dstar_collision"] = collision

        if right_pairs < min_right_pairs:
            min_right_pairs = right_pairs
            min_record = row_summary
        if dstar_right_pairs < min_dstar_right_pairs:
            min_dstar_right_pairs = dstar_right_pairs
            dstar_min_record = row_summary
        if min_edge_count is None or total_edges < min_edge_count:
            min_edge_count = total_edges
        if should_capture_thin_row and len(thin_records) < max_records:
            thin_records.append(row_summary)

        states = next_states
        if not states:
            break

        if progress_path and (row_index == 1 or row_index % 10000 == 0 or row_index == total_rows):
            write_json(
                progress_path,
                {
                    "status": "running",
                    "generated_at": datetime.now().isoformat(timespec="seconds"),
                    "p": p,
                    "q": q,
                    "rows_completed": row_index,
                    "rows_total": total_rows,
                    "percent": round(100.0 * row_index / total_rows, 2),
                    "min_active_B_after": min_active_b,
                    "min_right_pairs_after": min_right_pairs,
                    "min_Dstar_active_B_after": min_dstar_active_b,
                    "min_Dstar_pairs_after": min_dstar_right_pairs,
                    "min_transition_edge_count": min_edge_count,
                    "min_Dstar_phi_exact_color_m": min_dstar_phi_exact_color_m,
                },
            )

        min_active_b = min(min_active_b, active_b)
        min_dstar_active_b = min(min_dstar_active_b, dstar_active_b)

    return {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "p": p,
        "q": q,
        "rows_total": total_rows,
        "rows_completed": row_index if "row_index" in locals() else 0,
        "B_range": [b_min, b_max],
        "good_class_B": good_bs,
        "Dstar": sorted(dstar),
        "success": bool(states),
        "min_active_B_after": min_active_b,
        "min_right_pairs_after": min_right_pairs,
        "min_Dstar_active_B_after": min_dstar_active_b,
        "min_Dstar_pairs_after": min_dstar_right_pairs,
        "min_transition_edge_count": min_edge_count,
        "min_Dstar_phi_exact_color_m": min_dstar_phi_exact_color_m,
        "active_B_hist": {str(k): v for k, v in sorted(active_b_hist.items())},
        "right_pair_hist": {str(k): v for k, v in sorted(right_pair_hist.items())},
        "Dstar_pair_hist": {str(k): v for k, v in sorted(dstar_pair_hist.items())},
        "min_record": min_record,
        "Dstar_min_record": dstar_min_record,
        "Dstar_min_phi_record": min_dstar_phi_record,
        "thin_records": thin_records,
        "claim_boundary": [
            "one-edge deterministic layered transition diagnostic",
            "not a global proof",
            "not a replacement for the missing dynamic Hall-core lemma",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--output-json", required=True)
    parser.add_argument("--progress-json")
    parser.add_argument("--b-min", type=int, default=44)
    parser.add_argument("--b-max", type=int, default=74)
    parser.add_argument("--alpha", type=int, default=86)
    parser.add_argument("--radius", type=int, default=50)
    parser.add_argument("--width", type=int, default=30)
    parser.add_argument("--force-composite-xc", action="store_true")
    parser.add_argument("--skip-exact-hall", action="store_true", help="write thin rows only; compute exact Hall in postprocess")
    parser.add_argument("--dstar", default="-26,-22,-20,-10,-8,4,8,10,14")
    parser.add_argument("--max-records", type=int, default=20)
    args = parser.parse_args()

    dstar = {int(item) for item in args.dstar.split(",") if item.strip()}
    result = audit_manifest(
        read_json(Path(args.manifest)),
        b_min=args.b_min,
        b_max=args.b_max,
        alpha=args.alpha,
        radius=args.radius,
        width=args.width,
        force_composite=args.force_composite_xc,
        dstar=dstar,
        max_records=args.max_records,
        progress_path=Path(args.progress_json) if args.progress_json else None,
        skip_exact_hall=args.skip_exact_hall,
    )
    write_json(Path(args.output_json), result)
    print(
        json.dumps(
            {
                "output_json": args.output_json,
                "success": result["success"],
                "rows_completed": result["rows_completed"],
                "min_active_B_after": result["min_active_B_after"],
                "min_right_pairs_after": result["min_right_pairs_after"],
                "min_Dstar_active_B_after": result["min_Dstar_active_B_after"],
                "min_Dstar_pairs_after": result["min_Dstar_pairs_after"],
                "min_transition_edge_count": result["min_transition_edge_count"],
                "min_Dstar_phi_exact_color_m": result["min_Dstar_phi_exact_color_m"],
            },
            indent=2,
        )
    )
    return 0 if result["success"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
