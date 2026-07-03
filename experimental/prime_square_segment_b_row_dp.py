#!/usr/bin/env python3
"""
Row-DP monotone strip certificate for Segment B.

For each row y:
  1. Start with offsets reachable at row y.
  2. Close them under horizontal movement through Good vertices in that row.
  3. Lift to row y+1 at the same offsets that are Good on row y+1.

This is equivalent to monotone strip reachability, but exposes a row-by-row invariant
that is closer to a possible proof than unrestricted BFS.
"""

from __future__ import annotations

import argparse
import json
import math
from functools import lru_cache
from pathlib import Path

from sympy import isprime


@lru_cache(maxsize=None)
def prime(n: int) -> bool:
    return bool(isprime(n))


def good(x: int, y: int) -> bool:
    return math.gcd(x, y) == 1 and not (prime(x) and prime(y))


def components(good_offsets: set[int], min_d: int, max_d: int) -> list[tuple[int, int]]:
    comps = []
    d = min_d
    while d <= max_d:
        if d not in good_offsets:
            d += 1
            continue
        start = d
        while d + 1 <= max_d and d + 1 in good_offsets:
            d += 1
        comps.append((start, d))
        d += 1
    return comps


def horizontal_closure(reachable: set[int], good_offsets: set[int], min_d: int, max_d: int) -> set[int]:
    closed: set[int] = set()
    for a, b in components(good_offsets, min_d, max_d):
        if any(a <= r <= b for r in reachable):
            closed.update(range(a, b + 1))
    return closed


def run_row_dp(x_c: int, y0: int, y1: int, width: int) -> dict:
    min_d, max_d = -width, width
    reachable = {d for d in range(min_d, max_d + 1) if good(x_c + d, y0)}
    if not reachable:
        return {"success": False, "failure": "no_good_start", "failure_y": y0}

    stats = {
        "min_reachable_after_closure": None,
        "min_liftable": None,
        "max_components_per_row": 0,
        "min_largest_component": None,
        "max_abs_offset_reachable": 0,
        "rows_with_single_reachable": 0,
        "sample_rows": [],
    }

    for y in range(y0, y1):
        good_y = {d for d in range(min_d, max_d + 1) if good(x_c + d, y)}
        comps = components(good_y, min_d, max_d)
        closed = horizontal_closure(reachable, good_y, min_d, max_d)
        if not closed:
            return {
                "success": False,
                "failure": "no_horizontal_closure",
                "failure_y": y,
                "reachable_before": sorted(reachable),
                "good_offsets": sorted(good_y),
            }

        good_next = {d for d in range(min_d, max_d + 1) if good(x_c + d, y + 1)}
        liftable = closed & good_next
        if not liftable:
            return {
                "success": False,
                "failure": "no_liftable_offset",
                "failure_y": y,
                "closed_offsets": sorted(closed),
                "good_next_offsets": sorted(good_next),
            }

        largest = max((b - a + 1 for a, b in comps), default=0)
        stats["max_components_per_row"] = max(stats["max_components_per_row"], len(comps))
        stats["max_abs_offset_reachable"] = max(
            stats["max_abs_offset_reachable"],
            max(abs(d) for d in closed),
        )
        stats["min_reachable_after_closure"] = (
            len(closed)
            if stats["min_reachable_after_closure"] is None
            else min(stats["min_reachable_after_closure"], len(closed))
        )
        stats["min_liftable"] = (
            len(liftable)
            if stats["min_liftable"] is None
            else min(stats["min_liftable"], len(liftable))
        )
        stats["min_largest_component"] = (
            largest
            if stats["min_largest_component"] is None
            else min(stats["min_largest_component"], largest)
        )
        if len(closed) == 1:
            stats["rows_with_single_reachable"] += 1
        if len(stats["sample_rows"]) < 20:
            stats["sample_rows"].append({
                "y": y,
                "reachable_before": sorted(reachable),
                "closed_count": len(closed),
                "liftable_count": len(liftable),
                "components": comps,
            })

        reachable = liftable

    return {
        "success": True,
        "rows_tested": y1 - y0,
        "final_reachable": sorted(reachable),
        **stats,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cert-json", default="docs/front507_lemma3_detour_path_cert_2026-06-23.json")
    parser.add_argument("--pair-index", type=int, default=1)
    parser.add_argument("--width", type=int, default=30)
    parser.add_argument("--row-sample", type=int, default=0, help="0 means full segment")
    parser.add_argument("--output-json", default="docs/segment_b_row_dp_2026-06-24.json")
    args = parser.parse_args()

    cert = json.loads(Path(args.cert_json).read_text(encoding="utf-8"))
    pair = cert["results"][args.pair_index]
    p = int(pair["p"])
    p_next = int(pair["p_next"])
    x_c = int(pair["x_c"])
    y0 = p * p
    y1_full = p_next * p_next
    y1 = y1_full if args.row_sample == 0 else min(y0 + args.row_sample, y1_full)

    result = run_row_dp(x_c, y0, y1, args.width)
    result.update({
        "p": p,
        "p_next": p_next,
        "x_c": x_c,
        "y0": y0,
        "y1_full": y1_full,
        "y1_tested": y1,
        "rows_tested": y1 - y0,
        "width": args.width,
        "note": "Row-DP monotone strip certificate; empirical for selected segment.",
    })

    out = Path(args.output_json)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps({
        "output_json": args.output_json,
        "success": result["success"],
        "failure": result.get("failure"),
        "failure_y": result.get("failure_y"),
        "rows_tested": result["rows_tested"],
        "min_reachable_after_closure": result.get("min_reachable_after_closure"),
        "min_liftable": result.get("min_liftable"),
        "max_abs_offset_reachable": result.get("max_abs_offset_reachable"),
    }, indent=2))
    return 0 if result["success"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
