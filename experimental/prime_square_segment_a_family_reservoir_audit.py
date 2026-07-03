#!/usr/bin/env python3
"""Audit the good-class family reservoir R_y across one Segment A edge.

This targets the full-proof question: do all good-class B candidates ever get
simultaneously close to extinction during the row-DP interval?
"""

from __future__ import annotations

import argparse
import json
import math
from collections import Counter
from datetime import datetime
from pathlib import Path
from typing import Any

from sympy import factorint, isprime

from prime_square_segment_b_generic_sentinel_verify import DEFAULT_D
from prime_square_segment_b_row_dp import good, horizontal_closure


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def x_c_for(p: int, b: int, alpha: int, radius: int, force_composite: bool) -> tuple[int, int]:
    raw = alpha * p * p + b * p + radius + 1
    x_c = raw
    if force_composite:
        while isprime(x_c):
            x_c += 1
    return raw, x_c


def blocker_summary(x: int, y: int) -> dict[str, Any]:
    g = math.gcd(x, y)
    return {
        "gcd": g,
        "gcd_factors": {str(k): int(v) for k, v in factorint(g).items()} if g > 1 else {},
        "prime_prime": bool(isprime(x) and isprime(y)),
    }


def blocker_key(x: int, y: int) -> str:
    g = math.gcd(x, y)
    if g > 1:
        primes = sorted(int(p) for p in factorint(g))
        return "gcd:" + ",".join(str(p) for p in primes)
    if isprime(x) and isprime(y):
        return "prime_prime"
    return "survive"


def blocker_keys(x: int, y: int) -> set[str]:
    keys: set[str] = set()
    g = math.gcd(x, y)
    if g > 1:
        for prime in sorted(int(p) for p in factorint(g)):
            keys.add(f"gcd:{prime}")
    if isprime(x) and isprime(y):
        keys.add("prime_prime")
    return keys


def collision_lower_bound(vertices: list[tuple[int, int, int]], n: int) -> dict[str, Any]:
    blockers = {idx: blocker_keys(x, n) for idx, _b, _d, x in vertices}
    degrees = {idx: 0 for idx, *_ in vertices}
    blocker_sizes: Counter[str] = Counter()
    for keys in blockers.values():
        for key in keys:
            blocker_sizes[key] += 1
    for idx, *_ in vertices:
        neighbors: set[int] = set()
        for key in blockers[idx]:
            for other, other_keys in blockers.items():
                if other != idx and key in other_keys:
                    neighbors.add(other)
        degrees[idx] = len(neighbors)
    caro_wei = sum(1.0 / (degrees[idx] + 1) for idx in degrees)
    omega_n = len(factorint(n))
    m_n = omega_n + (1 if isprime(n) else 0)
    return {
        "n": n,
        "vertex_count": len(vertices),
        "omega_n": omega_n,
        "n_is_prime": bool(isprime(n)),
        "m_n": m_n,
        "caro_wei_alpha_lower_bound": caro_wei,
        "phi_caro_wei_lower_bound": caro_wei - m_n,
        "max_collision_degree": max(degrees.values()) if degrees else 0,
        "degree_hist": {str(k): v for k, v in sorted(Counter(degrees.values()).items())},
        "blocker_size_hist": {str(k): v for k, v in sorted(Counter(blocker_sizes.values()).items())},
        "top_blockers": [
            {"blocker": key, "size": size}
            for key, size in sorted(blocker_sizes.items(), key=lambda item: (-item[1], item[0]))[:20]
        ],
    }


def exact_independence_number(vertices: list[tuple[int, int, int, int]], n: int) -> int | None:
    if len(vertices) > 40:
        return None
    blockers = {idx: blocker_keys(x, n) for idx, _b, _d, x in vertices}
    adjacency = {idx: set() for idx, *_ in vertices}
    ids = [idx for idx, *_ in vertices]
    for i, left in enumerate(ids):
        for right in ids[i + 1 :]:
            if blockers[left] & blockers[right]:
                adjacency[left].add(right)
                adjacency[right].add(left)

    best = 0

    def bronk(candidates: set[int], chosen_size: int) -> None:
        nonlocal best
        if not candidates:
            best = max(best, chosen_size)
            return
        if chosen_size + len(candidates) <= best:
            return
        vertex = min(candidates, key=lambda item: len(adjacency[item] & candidates))
        bronk(candidates - {vertex}, chosen_size)
        bronk(candidates - adjacency[vertex] - {vertex}, chosen_size + 1)

    bronk(set(ids), 0)
    return best


def audit_edge(
    p: int,
    q: int,
    b_min: int,
    b_max: int,
    alpha: int,
    radius: int,
    width: int,
    force_composite: bool,
    event_threshold_active_b: int,
    event_threshold_pairs: int,
    max_events: int,
    progress_path: Path | None,
    skip_phi_scan: bool,
    skip_min_collision: bool,
) -> dict[str, Any]:
    bank_set = set(DEFAULT_D)
    y_start = p * p
    y_end = q * q
    good_bs = [b for b in range(b_min, b_max + 1) if b % 2 == 0 and b % 3 == p % 3]

    states: dict[int, set[int]] = {}
    shifted_bs: list[dict[str, int]] = []
    meta: dict[int, dict[str, int]] = {}
    for b in good_bs:
        raw, x_c = x_c_for(p, b, alpha, radius, force_composite)
        start = {d for d in bank_set if good(x_c + d, y_start)}
        if start:
            states[b] = start
        if x_c != raw:
            shifted_bs.append({"B": b, "raw_x_c": raw, "x_c": x_c, "force_shift": x_c - raw})
        meta[b] = {"raw_x_c": raw, "x_c": x_c, "force_shift": x_c - raw}

    active_hist: Counter[int] = Counter()
    pair_hist: Counter[int] = Counter()
    initial_nonempty_b = len(states)
    min_active_b = len(states)
    min_total_pairs = sum(len(v) for v in states.values())
    min_lift_pairs = min_total_pairs
    initial_classes = {(b % 30, d % 30) for b, bank in states.items() for d in bank}
    min_survivor_classes_mod30 = len(initial_classes)
    min_phi_caro_wei = None
    min_phi_alpha_exact = None
    min_phi_record: dict[str, Any] | None = None
    first_zero_y = None
    events: list[dict[str, Any]] = []
    min_record: dict[str, Any] | None = None

    total_rows = y_end - y_start
    for row_index, y in enumerate(range(y_start, y_end), start=1):
        next_states: dict[int, set[int]] = {}
        row_details = []
        total_lift_pairs = 0

        for b, bank in states.items():
            x_c = meta[b]["x_c"]
            good_y = {d for d in range(-width, width + 1) if good(x_c + d, y)}
            closure = horizontal_closure(bank, good_y, -width, width)
            good_next_bank = {d for d in bank_set if good(x_c + d, y + 1)}
            next_bank = closure & good_next_bank
            total_lift_pairs += len(next_bank)
            if next_bank:
                next_states[b] = next_bank
            row_details.append(
                {
                    "B": b,
                    "bank_before": sorted(bank),
                    "bank_before_size": len(bank),
                    "closure": sorted(closure),
                    "closure_size": len(closure),
                    "good_next_bank": sorted(good_next_bank),
                    "good_next_bank_size": len(good_next_bank),
                    "next_bank": sorted(next_bank),
                    "next_bank_size": len(next_bank),
                }
            )

        active_b = len(next_states)
        total_pairs = sum(len(v) for v in next_states.values())
        survivor_classes_mod30 = {(b % 30, d % 30) for b, bank in next_states.items() for d in bank}
        survivor_class_count_mod30 = len(survivor_classes_mod30)
        active_hist[active_b] += 1
        pair_hist[total_pairs] += 1
        min_active_b = min(min_active_b, active_b)
        min_total_pairs = min(min_total_pairs, total_pairs)
        min_lift_pairs = min(min_lift_pairs, total_lift_pairs)
        min_survivor_classes_mod30 = min(min_survivor_classes_mod30, survivor_class_count_mod30)

        if active_b == 0 and first_zero_y is None:
            first_zero_y = y

        if min_record is None or total_pairs < min_record["total_pairs_after"]:
            min_collision = None
            if not skip_min_collision:
                min_lift_vertices: list[tuple[int, int, int, int]] = []
                for detail in row_details:
                    b = detail["B"]
                    x_c = meta[b]["x_c"]
                    for d in set(detail["closure"]) & bank_set:
                        min_lift_vertices.append((len(min_lift_vertices), b, d, x_c + d))
                min_collision = collision_lower_bound(min_lift_vertices, y + 1)
                alpha_exact = exact_independence_number(min_lift_vertices, y + 1)
                if alpha_exact is not None:
                    min_collision["alpha_exact"] = alpha_exact
                    min_collision["phi_alpha_exact"] = alpha_exact - min_collision["m_n"]
            min_record = {
                "y": y,
                "next_y": y + 1,
                "active_b_after": active_b,
                "total_pairs_after": total_pairs,
                "survivor_class_count_mod30": survivor_class_count_mod30,
                "survivor_classes_mod30": sorted([list(item) for item in survivor_classes_mod30]),
                "surviving_B": sorted(next_states),
                "row_details": row_details,
            }
            if min_collision is not None:
                min_record["collision_lower_bound"] = min_collision

        if not skip_phi_scan and (len(events) < max_events or active_b <= event_threshold_active_b or total_pairs <= event_threshold_pairs):
            lift_vertices_for_phi: list[tuple[int, int, int, int]] = []
            for detail in row_details:
                b = detail["B"]
                x_c = meta[b]["x_c"]
                for d in set(detail["closure"]) & bank_set:
                    lift_vertices_for_phi.append((len(lift_vertices_for_phi), b, d, x_c + d))
            phi_collision = collision_lower_bound(lift_vertices_for_phi, y + 1)
            alpha_exact = exact_independence_number(lift_vertices_for_phi, y + 1)
            if alpha_exact is not None:
                phi_collision["alpha_exact"] = alpha_exact
                phi_collision["phi_alpha_exact"] = alpha_exact - phi_collision["m_n"]
                if min_phi_alpha_exact is None or phi_collision["phi_alpha_exact"] < min_phi_alpha_exact:
                    min_phi_alpha_exact = phi_collision["phi_alpha_exact"]
            if min_phi_caro_wei is None or phi_collision["phi_caro_wei_lower_bound"] < min_phi_caro_wei:
                min_phi_caro_wei = phi_collision["phi_caro_wei_lower_bound"]
                min_phi_record = {
                    "y": y,
                    "next_y": y + 1,
                    "active_b_after": active_b,
                    "total_pairs_after": total_pairs,
                    "collision_lower_bound": phi_collision,
                }

        if (
            (active_b <= event_threshold_active_b or total_pairs <= event_threshold_pairs)
            and len(events) < max_events
        ):
            blocker_hist: Counter[str] = Counter()
            lift_vertices: list[tuple[int, int, int, int]] = []
            for detail in row_details:
                b = detail["B"]
                x_c = meta[b]["x_c"]
                for d in set(detail["closure"]) & bank_set:
                    lift_vertices.append((len(lift_vertices), b, d, x_c + d))
                    if d not in set(detail["good_next_bank"]):
                        bs = blocker_summary(x_c + d, y + 1)
                        for prime, exp in bs["gcd_factors"].items():
                            blocker_hist[prime] += exp
                        if bs["prime_prime"]:
                            blocker_hist["prime_prime"] += 1
            events.append(
                {
                    "y": y,
                    "next_y": y + 1,
                    "active_b_after": active_b,
                    "total_pairs_after": total_pairs,
                    "survivor_class_count_mod30": survivor_class_count_mod30,
                    "survivor_classes_mod30": sorted([list(item) for item in survivor_classes_mod30]),
                    "total_lift_pairs_before_drop": total_lift_pairs,
                    "surviving_B": sorted(next_states),
                    "blocker_factor_hist": {str(k): v for k, v in sorted(blocker_hist.items(), key=lambda item: str(item[0]))},
                    "collision_lower_bound": collision_lower_bound(lift_vertices, y + 1),
                    "row_details": row_details,
                }
            )

        states = next_states
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
                    "current_y": y,
                    "min_active_B_after_step": min_active_b,
                    "min_total_pairs_after_step": min_total_pairs,
                    "min_survivor_classes_mod30": min_survivor_classes_mod30,
                    "min_phi_caro_wei_lower_bound": min_phi_caro_wei,
                    "min_phi_alpha_exact": min_phi_alpha_exact,
                },
            )
        if not states:
            break

    return {
        "p": p,
        "q": q,
        "gap": q - p,
        "B_range": [b_min, b_max],
        "good_B": good_bs,
        "initial_active_B": len(good_bs),
        "initial_nonempty_B": initial_nonempty_b,
        "force_composite_xc": force_composite,
        "shifted_good_B_count": len(shifted_bs),
        "shifted_good_B": shifted_bs,
        "success": first_zero_y is None,
        "first_zero_y": first_zero_y,
        "rows_tested": (first_zero_y - y_start + 1) if first_zero_y is not None else y_end - y_start,
        "min_active_B_after_step": min_active_b,
        "min_total_pairs_after_step": min_total_pairs,
        "min_lift_pairs_before_drop": min_lift_pairs,
        "min_survivor_classes_mod30": min_survivor_classes_mod30,
        "min_phi_caro_wei_lower_bound": min_phi_caro_wei,
        "min_phi_alpha_exact": min_phi_alpha_exact,
        "min_phi_record": min_phi_record,
        "active_B_hist": {str(k): v for k, v in sorted(active_hist.items())},
        "total_pairs_hist": {str(k): v for k, v in sorted(pair_hist.items())},
        "min_record": min_record,
        "final_states": {str(k): sorted(v) for k, v in sorted(states.items())},
        "events": events,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--output-json", required=True)
    ap.add_argument("--progress-json", default="")
    ap.add_argument("--b-min", type=int, default=44)
    ap.add_argument("--b-max", type=int, default=74)
    ap.add_argument("--alpha", type=int, default=86)
    ap.add_argument("--radius", type=int, default=50)
    ap.add_argument("--width", type=int, default=30)
    force_group = ap.add_mutually_exclusive_group()
    force_group.add_argument("--force-composite-xc", dest="force_composite_xc", action="store_true")
    force_group.add_argument("--raw-xc", dest="force_composite_xc", action="store_false")
    ap.set_defaults(force_composite_xc=True)
    ap.add_argument("--event-threshold-active-b", type=int, default=2)
    ap.add_argument("--event-threshold-pairs", type=int, default=8)
    ap.add_argument("--max-events", type=int, default=30)
    ap.add_argument("--skip-phi-scan", action="store_true")
    ap.add_argument("--skip-min-collision", action="store_true")
    args = ap.parse_args()

    manifest_path = Path(args.manifest)
    manifest = read_json(manifest_path)
    payload = audit_edge(
        p=int(manifest["p"]),
        q=int(manifest["q"]),
        b_min=args.b_min,
        b_max=args.b_max,
        alpha=args.alpha,
        radius=args.radius,
        width=args.width,
        force_composite=args.force_composite_xc,
        event_threshold_active_b=args.event_threshold_active_b,
        event_threshold_pairs=args.event_threshold_pairs,
        max_events=args.max_events,
        progress_path=Path(args.progress_json) if args.progress_json else None,
        skip_phi_scan=args.skip_phi_scan,
        skip_min_collision=args.skip_min_collision,
    )
    payload.update(
        {
            "generated_at": datetime.now().isoformat(timespec="seconds"),
            "claim_boundary": [
                "family reservoir diagnostic only",
                "not a global Segment A proof",
            ],
            "manifest_path": str(manifest_path).replace("\\", "/"),
        }
    )
    write_json(Path(args.output_json), payload)
    print(
        json.dumps(
            {
                "output_json": args.output_json,
                "success": payload["success"],
                "min_active_B_after_step": payload["min_active_B_after_step"],
                "min_total_pairs_after_step": payload["min_total_pairs_after_step"],
                "min_survivor_classes_mod30": payload["min_survivor_classes_mod30"],
                "min_phi_caro_wei_lower_bound": payload["min_phi_caro_wei_lower_bound"],
                "min_phi_alpha_exact": payload["min_phi_alpha_exact"],
                "shifted_good_B_count": payload["shifted_good_B_count"],
                "event_count": len(payload["events"]),
            },
            indent=2,
        )
    )
    return 0 if payload["success"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
