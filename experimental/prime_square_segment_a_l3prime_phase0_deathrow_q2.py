#!/usr/bin/env python3
"""L3' Phase-0 non-vacuous Q2 on death-bearing rows.

This runner measures reachable-bank lane shape exactly at lane-death events.
It is intentionally single-process and small-scope: the point is to avoid the
previous vacuous wide-window Q2 by working in a regime where deaths occur.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import time
from datetime import datetime
from pathlib import Path
from typing import Any

from sympy import factorint

from prime_square_segment_a_layered_transition_audit import component_id_map, x_c_for
from prime_square_segment_a_phase0_audit import shape_stats
from prime_square_segment_b_generic_sentinel_verify import DEFAULT_D
from prime_square_segment_b_row_dp import prime


D_STAR = {-26, -22, -20, -10, -8, 4, 8, 10, 14}


def parse_bank(text: str) -> list[int]:
    return [int(part.strip()) for part in text.split(",") if part.strip()]


def now_iso() -> str:
    return datetime.now().isoformat(timespec="seconds")


def atomic_write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    with tmp.open("w", encoding="utf-8", newline="\n") as f:
        f.write(text)
        f.flush()
        os.fsync(f.fileno())
    os.replace(tmp, path)


def write_json(path: Path, payload: dict[str, Any]) -> None:
    atomic_write_text(path, json.dumps(payload, indent=2, sort_keys=True) + "\n")


def initial_good_bs(p: int, b_min: int, w_lanes: int) -> list[int]:
    b = b_min
    while b % 2 != 0 or b % 3 != p % 3:
        b += 1
    return [b + 6 * i for i in range(w_lanes)]


def encode_states(states: dict[int, set[int]]) -> dict[str, list[int]]:
    return {str(b): sorted(v) for b, v in sorted(states.items())}


def decode_states(raw: dict[str, list[int]]) -> dict[int, set[int]]:
    return {int(b): {int(d) for d in values} for b, values in raw.items()}


def good_offsets_for_y(
    *,
    x_c: int,
    y: int,
    offsets: range,
    x_prime_by_d: dict[int, bool],
    y_is_prime: bool,
) -> set[int]:
    if y_is_prime:
        return {d for d in offsets if math.gcd(x_c + d, y) == 1 and not x_prime_by_d[d]}
    return {d for d in offsets if math.gcd(x_c + d, y) == 1}


def transition_fast(
    *,
    states: dict[int, set[int]],
    x_meta: dict[int, int],
    x_prime: dict[int, dict[int, bool]],
    y: int,
    width: int,
    bank_set: set[int],
) -> tuple[dict[int, set[int]], dict[int, dict[str, Any]]]:
    next_states: dict[int, set[int]] = {}
    deaths: dict[int, dict[str, Any]] = {}
    full_offsets = range(-width, width + 1)
    y_is_prime = prime(y)
    y_next_is_prime = prime(y + 1)
    for b, bank in states.items():
        x_c = x_meta[b]
        good_y = good_offsets_for_y(
            x_c=x_c,
            y=y,
            offsets=full_offsets,
            x_prime_by_d=x_prime[b],
            y_is_prime=y_is_prime,
        )
        comp = component_id_map(good_y, -width, width)
        left_components = {comp[d] for d in bank if d in comp}
        good_next = good_offsets_for_y(
            x_c=x_c,
            y=y + 1,
            offsets=bank_set,
            x_prime_by_d=x_prime[b],
            y_is_prime=y_next_is_prime,
        )
        next_bank = {e for e in good_next if e in comp and comp[e] in left_components}
        if next_bank:
            next_states[b] = next_bank
        else:
            deaths[b] = {
                "B": b,
                "death_y": y + 1,
                "transition_y": y,
                "n": y + 1,
                "n_factors": sorted(int(p0) for p0 in factorint(y + 1)),
                "bank_before": sorted(bank),
                "Dstar_before_count": len(bank & D_STAR),
                "good_y_count": len(good_y),
                "good_next_count": len(good_next),
                "left_component_count": len(left_components),
            }
    return next_states, deaths


def lane_snapshot(good_bs: list[int], states: dict[int, set[int]]) -> dict[str, Any]:
    live_bs = [b for b in good_bs if b in states]
    dstar_bs = [b for b in live_bs if states[b] & D_STAR]
    dead_bs = [b for b in good_bs if b not in states]
    dstar_counts = {str(b): len(states[b] & D_STAR) for b in live_bs}
    bank_sizes = {str(b): len(states[b]) for b in live_bs}
    live_shape = shape_stats(good_bs, live_bs)
    dstar_shape = shape_stats(good_bs, dstar_bs)
    return {
        "live_B": live_bs,
        "dead_B": dead_bs,
        "dstar_nonempty_B": dstar_bs,
        "live_shape": live_shape,
        "dstar_shape": dstar_shape,
        "living_Dstar_pairs": sum(int(v) for v in dstar_counts.values()),
        "living_Dstar_nonempty_B_count": len(dstar_bs),
        "Dstar_counts": dstar_counts,
        "bank_sizes": bank_sizes,
    }


def update_shape_extrema(extrema: dict[str, Any], y: int, y_start: int, states: dict[int, set[int]], good_bs: list[int]) -> None:
    snap = lane_snapshot(good_bs, states)
    for key, shape_key, better in [
        ("min_live_shape", "live_shape", "min"),
        ("min_dstar_shape", "dstar_shape", "min"),
        ("max_live_fragments", "live_shape", "max_frag"),
        ("max_dstar_fragments", "dstar_shape", "max_frag"),
    ]:
        stats = snap[shape_key]
        rec = {"y": y, "row_offset": y - y_start, **stats}
        old = extrema.get(key)
        if old is None:
            extrema[key] = rec
        elif better == "min" and (
            stats["maxrun_index"] < old["maxrun_index"]
            or (stats["maxrun_index"] == old["maxrun_index"] and stats["gapratio_index"] < old["gapratio_index"])
        ):
            extrema[key] = rec
        elif better == "max_frag" and (
            stats["fragment_count_index"] > old["fragment_count_index"]
            or (
                stats["fragment_count_index"] == old["fragment_count_index"]
                and stats["gapratio_index"] < old["gapratio_index"]
            )
        ):
            extrema[key] = rec
    old_pairs = extrema.get("min_living_Dstar_pairs")
    if old_pairs is None or snap["living_Dstar_pairs"] < old_pairs:
        extrema["min_living_Dstar_pairs"] = snap["living_Dstar_pairs"]
        extrema["min_living_Dstar_pairs_record"] = {"y": y, "row_offset": y - y_start, **snap}


def compatible_progress(payload: dict[str, Any], args: argparse.Namespace, rows_total: int) -> bool:
    bank_values = sorted(parse_bank(args.bank))
    return (
        payload.get("runner") == "prime_square_segment_a_l3prime_phase0_deathrow_q2.v1"
        and payload.get("status") == "running"
        and payload.get("p") == args.p
        and payload.get("q") == args.q
        and payload.get("W_lanes") == args.w_lanes
        and payload.get("rows_total") == rows_total
        and payload.get("alpha") == args.alpha
        and payload.get("radius") == args.radius
        and payload.get("width") == args.width
        and payload.get("b_min") == args.b_min
        and payload.get("force_composite_xc") == bool(args.force_composite_xc)
        and payload.get("bank") == bank_values
        and "states" in payload
    )


def progress_payload(
    *,
    args: argparse.Namespace,
    rows_done: int,
    rows_total: int,
    y_start: int,
    states: dict[int, set[int]],
    good_bs: list[int],
    death_events: list[dict[str, Any]],
    extrema: dict[str, Any],
    started_at: float,
    status: str,
) -> dict[str, Any]:
    return {
        "runner": "prime_square_segment_a_l3prime_phase0_deathrow_q2.v1",
        "status": status,
        "p": args.p,
        "q": args.q,
        "W_lanes": args.w_lanes,
        "alpha": args.alpha,
        "radius": args.radius,
        "width": args.width,
        "b_min": args.b_min,
        "force_composite_xc": bool(args.force_composite_xc),
        "bank": sorted(parse_bank(args.bank)),
        "rows_done": rows_done,
        "rows_total": rows_total,
        "percent": round(100.0 * rows_done / rows_total, 6) if rows_total else 100.0,
        "wall_elapsed_s": round(time.monotonic() - started_at, 3),
        "last_update_iso": now_iso(),
        "current_y": y_start + rows_done,
        "death_event_count": len(death_events),
        "distinct_death_rows": len({event["death_y"] for event in death_events}),
        "live_lane_count": len(states),
        "states": encode_states(states),
        "death_events": death_events,
        "shape_extrema": extrema,
    }


def load_progress(path: Path, args: argparse.Namespace, rows_total: int) -> tuple[int, dict[int, set[int]] | None, list[dict[str, Any]], dict[str, Any]]:
    if not path.exists():
        return 0, None, [], {}
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return 0, None, [], {"resume_warning": "progress_json_unreadable_started_from_zero"}
    if payload.get("status") == "complete":
        return rows_total, decode_states(payload.get("states", {})), payload.get("death_events", []), payload.get("shape_extrema", {})
    if not compatible_progress(payload, args, rows_total):
        return 0, None, [], {"resume_warning": "incompatible_progress_json_started_from_zero"}
    return (
        int(payload.get("rows_done", 0)),
        decode_states(payload.get("states", {})),
        list(payload.get("death_events", [])),
        dict(payload.get("shape_extrema", {})),
    )


def run(args: argparse.Namespace) -> dict[str, Any]:
    y_start = args.p * args.p
    y_end_full = args.q * args.q
    y_end = y_end_full if args.max_rows is None else min(y_end_full, y_start + args.max_rows)
    rows_total = y_end - y_start
    good_bs = initial_good_bs(args.p, args.b_min, args.w_lanes)
    bank_values = parse_bank(args.bank)
    bank_set = set(bank_values)
    if max(abs(d) for d in bank_set) > args.width:
        raise SystemExit("--width must cover --bank")

    x_meta: dict[int, int] = {}
    x_prime: dict[int, dict[int, bool]] = {}
    for b in good_bs:
        _raw, x_c = x_c_for(args.p, b, args.alpha, args.radius, args.force_composite_xc)
        x_meta[b] = x_c
        x_prime[b] = {d: prime(x_c + d) for d in range(-args.width, args.width + 1)}

    started_at = time.monotonic()
    progress_path = Path(args.progress_json) if args.progress_json else None
    rows_done = 0
    states: dict[int, set[int]] | None = None
    death_events: list[dict[str, Any]] = []
    extrema: dict[str, Any] = {}
    resume_info: dict[str, Any] = {}
    if progress_path:
        rows_done, states, death_events, extrema = load_progress(progress_path, args, rows_total)
        if states is not None:
            resume_info = {"resume_status": "resumed", "resume_rows_done": rows_done}

    if states is None:
        y_start_is_prime = prime(y_start)
        states = {}
        for b in good_bs:
            start = good_offsets_for_y(
                x_c=x_meta[b],
                y=y_start,
                offsets=bank_set,
                x_prime_by_d=x_prime[b],
                y_is_prime=y_start_is_prime,
            )
            if start:
                states[b] = start
        update_shape_extrema(extrema, y_start, y_start, states, good_bs)

    last_heartbeat = 0.0
    last_rows_done = rows_done
    last_progress_time = time.monotonic()

    def heartbeat(status: str) -> None:
        payload = progress_payload(
            args=args,
            rows_done=rows_done,
            rows_total=rows_total,
            y_start=y_start,
            states=states or {},
            good_bs=good_bs,
            death_events=death_events,
            extrema=extrema,
            started_at=started_at,
            status=status,
        )
        if progress_path:
            write_json(progress_path, payload)
        print(
            f"L3Q2 W={args.w_lanes} {payload['percent']}% "
            f"rows={rows_done}/{rows_total} deaths={len(death_events)} live={len(states or {})}",
            flush=True,
        )

    heartbeat("running")
    while rows_done < rows_total and states:
        if args.stop_after_deaths is not None and len(death_events) >= args.stop_after_deaths:
            break
        y = y_start + rows_done
        before_states = {b: set(bank) for b, bank in states.items()}
        before_snap = lane_snapshot(good_bs, before_states)
        next_states, deaths = transition_fast(
            states=states,
            x_meta=x_meta,
            x_prime=x_prime,
            y=y,
            width=args.width,
            bank_set=bank_set,
        )
        rows_done += 1
        states = next_states
        update_shape_extrema(extrema, y + 1, y_start, states, good_bs)
        if deaths:
            after_snap = lane_snapshot(good_bs, states)
            for b, rec in sorted(deaths.items()):
                event = {
                    **rec,
                    "row_offset": (y + 1) - y_start,
                    "interior_death": b not in (good_bs[0], good_bs[-1]),
                    "before_shape": before_snap,
                    "after_shape": after_snap,
                }
                death_events.append(event)
            if args.stop_after_deaths is not None and len(death_events) >= args.stop_after_deaths:
                update_shape_extrema(extrema, y + 1, y_start, states, good_bs)

        now = time.monotonic()
        if rows_done != last_rows_done:
            last_rows_done = rows_done
            last_progress_time = now
        elif now - last_progress_time > args.no_progress_timeout_s:
            heartbeat("stalled")
            raise SystemExit(f"no row progress for {args.no_progress_timeout_s}s")
        if now - last_heartbeat >= args.heartbeat_seconds or rows_done % args.heartbeat_rows == 0:
            heartbeat("running")
            last_heartbeat = now

    distinct_death_rows = sorted({event["death_y"] for event in death_events})
    interior_deaths = [event for event in death_events if event["interior_death"]]
    max_after_live_fragments = max((event["after_shape"]["live_shape"]["fragment_count_index"] for event in death_events), default=0)
    max_after_dstar_fragments = max((event["after_shape"]["dstar_shape"]["fragment_count_index"] for event in death_events), default=0)
    max_after_live_gapratio_loss = min((event["after_shape"]["live_shape"]["gapratio_index"] for event in death_events), default=1.0)
    non_vacuity_pass = (
        len(distinct_death_rows) >= args.min_distinct_death_rows
        and len(interior_deaths) >= args.min_interior_deaths
        and max_after_live_fragments >= args.min_separators + 1
    )
    if not death_events:
        verdict = "VACUOUS"
    elif max_after_live_fragments <= 1 and max_after_dstar_fragments <= 1:
        verdict = "ALIVE_INTERVAL_LIKE"
    else:
        verdict = "DEAD_FRAGMENTED"

    payload = {
        "generated_at": now_iso(),
        "runner": "prime_square_segment_a_l3prime_phase0_deathrow_q2.v1",
        "status": "complete",
        "p": args.p,
        "q": args.q,
        "W_lanes": args.w_lanes,
        "B_range": [good_bs[0], good_bs[-1]],
        "good_class_B": good_bs,
        "alpha": args.alpha,
        "radius": args.radius,
        "width": args.width,
        "b_min": args.b_min,
        "force_composite_xc": bool(args.force_composite_xc),
        "bank": sorted(bank_set),
        "rows_done": rows_done,
        "rows_total": rows_total,
        "rows_total_full_interval": y_end_full - y_start,
        "scan_complete_interval": rows_done == y_end_full - y_start,
        "percent": 100.0 if rows_done == rows_total else round(100.0 * rows_done / rows_total, 6),
        "wall_elapsed_s": round(time.monotonic() - started_at, 3),
        "death_event_count": len(death_events),
        "distinct_death_rows": distinct_death_rows,
        "interior_death_count": len(interior_deaths),
        "max_after_live_fragments": max_after_live_fragments,
        "max_after_dstar_fragments": max_after_dstar_fragments,
        "min_after_live_gapratio_on_death_rows": max_after_live_gapratio_loss,
        "non_vacuity_gate": {
            "pass": non_vacuity_pass,
            "min_distinct_death_rows": args.min_distinct_death_rows,
            "min_interior_deaths": args.min_interior_deaths,
            "min_separators": args.min_separators,
            "observed_distinct_death_rows": len(distinct_death_rows),
            "observed_interior_deaths": len(interior_deaths),
            "observed_separators_from_live_fragments": max(0, max_after_live_fragments - 1),
        },
        "verdict": verdict,
        "shape_extrema": extrema,
        "death_events": death_events,
        "final_states": encode_states(states or {}),
        "resume_info": resume_info,
    }
    if progress_path:
        write_json(progress_path, {**payload, "last_update_iso": now_iso()})
    return payload


def write_md(path: Path, payload: dict[str, Any]) -> None:
    lines = [
        "# L3' Phase-0 Death-Row Q2",
        "",
        f"- generated_at: `{payload['generated_at']}`",
        f"- status: `{payload['status']}`",
        f"- p: `{payload['p']}`",
        f"- q: `{payload['q']}`",
        f"- W_lanes: `{payload['W_lanes']}`",
        f"- B_range: `{payload['B_range']}`",
        f"- rows_done: `{payload['rows_done']}` / `{payload['rows_total']}`",
        f"- death_event_count: `{payload['death_event_count']}`",
        f"- distinct_death_rows: `{len(payload['distinct_death_rows'])}`",
        f"- interior_death_count: `{payload['interior_death_count']}`",
        f"- max_after_live_fragments: `{payload['max_after_live_fragments']}`",
        f"- max_after_dstar_fragments: `{payload['max_after_dstar_fragments']}`",
        f"- non_vacuity_gate: `{payload['non_vacuity_gate']['pass']}`",
        f"- verdict: `{payload['verdict']}`",
        "",
        "## Non-Vacuity Gate",
        "",
        "```json",
        json.dumps(payload["non_vacuity_gate"], indent=2, sort_keys=True),
        "```",
        "",
        "## Death Events",
        "",
        "| B | death_y | row_offset | interior | after live count | after live MAXRUN | after live GAPRATIO | after live fragments | after D* MAXRUN | after D* fragments |",
        "|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|",
    ]
    for event in payload["death_events"]:
        live = event["after_shape"]["live_shape"]
        dstar = event["after_shape"]["dstar_shape"]
        lines.append(
            f"| {event['B']} | {event['death_y']} | {event['row_offset']} | {event['interior_death']} | "
            f"{live['count']} | {live['maxrun_index']} | {live['gapratio_index']:.6f} | "
            f"{live['fragment_count_index']} | {dstar['maxrun_index']} | {dstar['fragment_count_index']} |"
        )
    lines.extend(
        [
            "",
            "## Shape Extrema",
            "",
            "```json",
            json.dumps(payload["shape_extrema"], indent=2, sort_keys=True),
            "```",
            "",
            "## Interpretation",
            "",
        ]
    )
    if payload["verdict"] == "VACUOUS":
        lines.append("This run is vacuous: it did not observe lane deaths.")
    elif payload["verdict"] == "ALIVE_INTERVAL_LIKE":
        lines.append("Death rows remained interval-like under the configured criteria.")
    else:
        lines.append(
            "Death rows fragmented the surviving lane set. This rules out the tested "
            "edge052/W=11 index-interval invariant as stated; it does not rule out a "
            "fragmented general DP argument."
        )
    lines.append("")
    atomic_write_text(path, "\n".join(lines))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--p", type=int, required=True)
    ap.add_argument("--q", type=int, required=True)
    ap.add_argument("--w-lanes", type=int, required=True)
    ap.add_argument("--b-min", type=int, default=44)
    ap.add_argument("--alpha", type=int, default=86)
    ap.add_argument("--radius", type=int, default=80)
    ap.add_argument("--width", type=int, default=80)
    ap.add_argument("--bank", default=",".join(str(d) for d in DEFAULT_D))
    ap.add_argument("--force-composite-xc", action="store_true")
    ap.add_argument("--max-rows", type=int)
    ap.add_argument("--heartbeat-seconds", type=float, default=2.0)
    ap.add_argument("--heartbeat-rows", type=int, default=25000)
    ap.add_argument("--no-progress-timeout-s", type=float, default=60.0)
    ap.add_argument("--min-distinct-death-rows", type=int, default=3)
    ap.add_argument("--min-interior-deaths", type=int, default=3)
    ap.add_argument("--min-separators", type=int, default=2)
    ap.add_argument("--stop-after-deaths", type=int)
    ap.add_argument("--progress-json")
    ap.add_argument("--output-json", required=True)
    ap.add_argument("--output-md", required=True)
    args = ap.parse_args()

    payload = run(args)
    write_json(Path(args.output_json), payload)
    write_md(Path(args.output_md), payload)
    print(f"WROTE {args.output_json}", flush=True)
    print(f"WROTE {args.output_md}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
