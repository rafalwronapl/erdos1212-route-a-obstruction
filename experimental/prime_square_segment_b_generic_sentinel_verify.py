#!/usr/bin/env python3
"""Generic sentinel-bank row-DP verifier for one proposed connector pair."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from prime_square_segment_b_row_dp import good, horizontal_closure


DEFAULT_D = [-26, -22, -20, -10, -8, -4, -2, 2, 4, 8, 10, 14, 20, 22, 28]


def parse_bank(text: str) -> list[int]:
    return [int(x) for x in text.split(",") if x.strip()]


def run(x_c: int, y_start: int, y_end: int, bank: list[int], width: int,
        weak_threshold: int, max_weak_events: int) -> dict:
    bank_set = set(bank)
    b = {d for d in bank_set if good(x_c + d, y_start)}
    if not b:
        return {"success": False, "failure": "empty_start", "failure_y": y_start}

    min_bank = len(b)
    max_bank = len(b)
    weak_events = []

    for y in range(y_start, y_end):
        good_y = {d for d in range(-width, width + 1) if good(x_c + d, y)}
        closure = horizontal_closure(b, good_y, -width, width)
        next_b = closure & {d for d in bank_set if good(x_c + d, y + 1)}
        if not next_b:
            return {
                "success": False,
                "failure": "empty_next_bank",
                "failure_y": y,
                "bank_before": sorted(b),
                "closure_size": len(closure),
            }
        b = next_b
        min_bank = min(min_bank, len(b))
        max_bank = max(max_bank, len(b))
        if len(b) <= weak_threshold and len(weak_events) < max_weak_events:
            weak_events.append({"y": y + 1, "bank_size": len(b), "bank": sorted(b)})

    return {
        "success": True,
        "rows_tested": y_end - y_start,
        "y_start": y_start,
        "y_end": y_end,
        "final_bank": sorted(b),
        "min_bank": min_bank,
        "max_bank": max_bank,
        "weak_threshold": weak_threshold,
        "weak_event_count_recorded": len(weak_events),
        "weak_events": weak_events,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--p", type=int, required=True)
    ap.add_argument("--q", type=int, required=True)
    ap.add_argument("--x-c", type=int, required=True)
    ap.add_argument("--width", type=int, default=30)
    ap.add_argument("--bank", default=",".join(str(x) for x in DEFAULT_D))
    ap.add_argument("--weak-threshold", type=int, default=5)
    ap.add_argument("--max-weak-events", type=int, default=50)
    ap.add_argument("--output-json", required=True)
    args = ap.parse_args()

    bank = parse_bank(args.bank)
    result = run(
        x_c=args.x_c,
        y_start=args.p * args.p,
        y_end=args.q * args.q,
        bank=bank,
        width=args.width,
        weak_threshold=args.weak_threshold,
        max_weak_events=args.max_weak_events,
    )
    result.update({
        "p": args.p,
        "q": args.q,
        "x_c": args.x_c,
        "width": args.width,
        "bank_D": bank,
    })
    Path(args.output_json).write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps({
        "output_json": args.output_json,
        "success": result.get("success"),
        "rows_tested": result.get("rows_tested"),
        "min_bank": result.get("min_bank"),
        "max_bank": result.get("max_bank"),
        "weak_event_count_recorded": result.get("weak_event_count_recorded"),
        "failure": result.get("failure"),
        "failure_y": result.get("failure_y"),
    }, indent=2))
    return 0 if result.get("success") else 1


if __name__ == "__main__":
    raise SystemExit(main())
