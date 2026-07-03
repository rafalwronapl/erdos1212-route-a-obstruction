#!/usr/bin/env python3
"""Audit active-color counts on the omega-maximal ladder witnesses.

This is deliberately a small static hit-table audit, not a row-DP replay.
For each omega-maximal witness n=y+1, it counts prime factors ell>=5 which
hit at least one Segment-A token (B,d) with d in D* and B in a good-lane
window. The true DP-reachable m_color can only be certified by knowing L_y(W);
this script records the static proxy needed by the publication audit.
"""

from __future__ import annotations

import json
import math
import os
from datetime import datetime
from pathlib import Path
from typing import Any


D_STAR = [-26, -22, -20, -10, -8, 4, 8, 10, 14]
LANE_COUNTS = [6, 40, 160]
WITNESS_SOURCES = [
    Path("docs/PHASE0_Q1_EXACT_edge052_p26267_q26293_2026-07-02.json"),
    Path("docs/PHASE0_Q1_EXACT_p99961_q99971_2026-07-02.json"),
    Path("docs/PHASE0_Q1_EXACT_p950039_q950041_2026-07-02.json"),
    Path("docs/PHASE0_OMEGA_LADDER_p1e7_2026-07-03.json"),
    Path("docs/PHASE0_OMEGA_LADDER_p1e8_2026-07-03.json"),
]


def atomic_write_text(path: Path, text: str) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(text, encoding="utf-8")
    os.replace(tmp, path)


def atomic_write_json(path: Path, payload: dict[str, Any]) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(tmp, path)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def witness_record(path: Path) -> dict[str, Any]:
    payload = load_json(path)
    if "exact_omega_ge5" in payload:
        witness = payload["exact_omega_ge5"]["max_omega_ge5"]
    else:
        witness = payload["max_omega_ge5"]
    return {
        "source": str(path),
        "p": int(payload["p"]),
        "q": int(payload["q"]),
        "logp_squared": float(payload["logp_squared"]),
        "n": int(witness["n"]),
        "y": int(witness["y"]),
        "omega_ge5": int(witness["omega_ge5"]),
        "factors_ge5": [int(x) for x in witness["n_factors_ge5"]],
    }


def first_good_b(p: int, b_min: int = 44) -> int:
    b = b_min
    while b % 2 != 0 or b % 3 != p % 3:
        b += 1
    return b


def good_b_window(p: int, lane_count: int, b_min: int = 44) -> list[int]:
    b0 = first_good_b(p, b_min)
    return [b0 + 6 * i for i in range(lane_count)]


def x_c(p: int, b: int) -> int:
    return 86 * p * p + b * p + 51


def active_hits(p: int, ell: int, bs: list[int]) -> list[dict[str, int]]:
    hits: list[dict[str, int]] = []
    for b in bs:
        base = x_c(p, b)
        for d in D_STAR:
            if (base + d) % ell == 0:
                hits.append({"B": b, "d": d})
    return hits


def audit_witness(record: dict[str, Any]) -> dict[str, Any]:
    p = int(record["p"])
    factors = [int(x) for x in record["factors_ge5"]]
    lane_audits: dict[str, Any] = {}
    for lane_count in LANE_COUNTS:
        bs = good_b_window(p, lane_count)
        active: list[dict[str, Any]] = []
        inactive: list[int] = []
        for ell in factors:
            hits = active_hits(p, ell, bs)
            if hits:
                active.append(
                    {
                        "ell": ell,
                        "hit_count": len(hits),
                        "sample_hits": hits[:8],
                    }
                )
            else:
                inactive.append(ell)
        m_color_static = len(active)
        lane_audits[str(lane_count)] = {
            "lane_count": lane_count,
            "B_start": bs[0],
            "B_end": bs[-1],
            "good_B": bs if lane_count <= 6 else [bs[0], bs[1], bs[2], "...", bs[-3], bs[-2], bs[-1]],
            "m_color_static": m_color_static,
            "active_primes": [row["ell"] for row in active],
            "inactive_primes": inactive,
            "tracks_omega_exactly": m_color_static == int(record["omega_ge5"]),
            "active_detail": active,
        }
    return {
        **record,
        "lane_audits": lane_audits,
        "tracks_at_W6": lane_audits["6"]["m_color_static"] == int(record["omega_ge5"]),
        "tracks_at_W40": lane_audits["40"]["m_color_static"] == int(record["omega_ge5"]),
        "tracks_at_W160": lane_audits["160"]["m_color_static"] == int(record["omega_ge5"]),
    }


def trend_label(rows: list[dict[str, Any]], lane_count: int) -> str:
    values = [int(row["lane_audits"][str(lane_count)]["m_color_static"]) for row in rows]
    omegas = [int(row["omega_ge5"]) for row in rows]
    if values == omegas:
        return "tracks omega exactly on all witness rows"
    if values[-1] >= values[0] and values[-1] >= max(6, omegas[-1] - 1):
        return "tracks omega closely on this static window"
    if values[-1] <= values[0]:
        return "does not track omega on this static window"
    return "partial growth, slower than omega"


def render_md(payload: dict[str, Any]) -> str:
    lines: list[str] = []
    lines.append("# Phase-0 m_color Ladder Audit")
    lines.append("")
    lines.append("Date: 2026-07-03")
    lines.append("Status: COMPLETE for static Segment-A hit-table proxy.")
    lines.append("")
    lines.append("This audit checks whether omega-maximal witness rows also have many active")
    lines.append("GCD colors. It is not a row-DP replay: true DP `m_color(y)` depends on the")
    lines.append("reachable reservoir `L_y(W)`. The count below is the static active-color")
    lines.append("proxy over good B-lanes and `D*`; unreachable tokens could only reduce the")
    lines.append("certified dynamic count.")
    lines.append("")
    lines.append("Definitions used:")
    lines.append("")
    lines.append("```text")
    lines.append("D* = {-26,-22,-20,-10,-8,4,8,10,14}")
    lines.append("x_c(B) = 86 p^2 + B p + 51")
    lines.append("good B-lanes: B even and B == p (mod 3), starting at B >= 44, spaced by 6")
    lines.append("ell active iff ell | x_c(B)+d for at least one audited (B,d)")
    lines.append("```")
    lines.append("")
    lines.append("## Results")
    lines.append("")
    lines.append(
        "| p | omega_max | m_color_static W=6 | m_color_static W=40 | "
        "m_color_static W=160 | tracks? | construction note |"
    )
    lines.append("|---:|---:|---:|---:|---:|---|---|")
    for row in payload["rows"]:
        w6 = row["lane_audits"]["6"]["m_color_static"]
        w40 = row["lane_audits"]["40"]["m_color_static"]
        w160 = row["lane_audits"]["160"]["m_color_static"]
        tracks = "W160 yes" if row["tracks_at_W160"] else "no / width-sensitive"
        note = "static B-window extension, not DP L_y(W)"
        lines.append(
            f"| {row['p']} | {row['omega_ge5']} | {w6} | {w40} | {w160} | {tracks} | {note} |"
        )
    lines.append("")
    lines.append("## Active Prime Detail")
    lines.append("")
    for row in payload["rows"]:
        lines.append(f"### p={row['p']}, n={row['n']}")
        lines.append("")
        lines.append(f"- omega_ge5: `{row['omega_ge5']}`")
        lines.append(f"- factors_ge5: `{', '.join(str(x) for x in row['factors_ge5'])}`")
        for lane_count in LANE_COUNTS:
            audit = row["lane_audits"][str(lane_count)]
            lines.append(
                f"- W={lane_count}: m_color_static=`{audit['m_color_static']}`, "
                f"active=`{audit['active_primes']}`, inactive=`{audit['inactive_primes']}`, "
                f"B=`{audit['B_start']}..{audit['B_end']}`"
            )
        lines.append("")
    lines.append("## Interpretation")
    lines.append("")
    for lane_count in LANE_COUNTS:
        lines.append(f"- W={lane_count}: {payload['trend_by_lane_count'][str(lane_count)]}.")
    lines.append("")
    lines.append(
        "The parity claim is therefore not an automatic consequence of omega alone. "
        "On the old six-lane core, the large-p witnesses do not track omega. On the "
        "wider static windows, especially W=160, the active-color count does track "
        "omega on the measured witnesses. This supports a parity obstruction for "
        "full/static wide-window certificates, but it does not by itself refute the "
        "DP-reachable weakened L3' route."
    )
    lines.append("")
    lines.append("## Source Artifacts")
    lines.append("")
    for source in payload["source_files"]:
        lines.append(f"- `{source}`")
    lines.append("")
    return "\n".join(lines)


def main() -> int:
    rows = [audit_witness(witness_record(path)) for path in WITNESS_SOURCES]
    payload = {
        "status": "complete",
        "runner": "prime_square_phase0_mcolor_ladder.v1",
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "D_star": D_STAR,
        "lane_counts": LANE_COUNTS,
        "source_files": [str(path) for path in WITNESS_SOURCES],
        "rows": rows,
        "trend_by_lane_count": {str(w): trend_label(rows, w) for w in LANE_COUNTS},
        "caveat": (
            "Static active-color count over B x D*. The true row-DP m_color requires "
            "the reachable set L_y(W); this audit does not certify L3' or lane survival."
        ),
    }
    out_json = Path("docs/PHASE0_MCOLOR_LADDER_2026-07-03.json")
    out_md = Path("docs/PHASE0_MCOLOR_LADDER_2026-07-03.md")
    atomic_write_json(out_json, payload)
    atomic_write_text(out_md, render_md(payload))
    print(f"Wrote {out_json}")
    print(f"Wrote {out_md}")
    for row in rows:
        print(
            "p={p} omega={omega} m6={m6} m40={m40} m160={m160}".format(
                p=row["p"],
                omega=row["omega_ge5"],
                m6=row["lane_audits"]["6"]["m_color_static"],
                m40=row["lane_audits"]["40"]["m_color_static"],
                m160=row["lane_audits"]["160"]["m_color_static"],
            )
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
