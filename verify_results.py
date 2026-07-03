#!/usr/bin/env python3
"""Lightweight consistency checks for the public Route-A obstruction packet."""

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"


def load(name: str) -> dict:
    return json.loads((RESULTS / name).read_text(encoding="utf-8"))


def check_mcolor_ladder() -> None:
    data = load("PHASE0_MCOLOR_LADDER_2026-07-03.json")
    rows = data["rows"]
    expected_p = [26267, 99961, 950039, 10000019, 100000007]
    expected_omega = [7, 7, 8, 10, 11]
    expected_w40 = [7, 7, 8, 10, 11]
    expected_w160 = [7, 7, 8, 10, 11]

    assert [r["p"] for r in rows] == expected_p
    assert [r["omega_ge5"] for r in rows] == expected_omega
    assert [r["lane_audits"]["40"]["m_color_static"] for r in rows] == expected_w40
    assert [r["lane_audits"]["160"]["m_color_static"] for r in rows] == expected_w160


def check_l3_deathrow() -> None:
    data = load("L3PRIME_PHASE0_DEATHROW_Q2_edge052_W11_2026-07-03.json")
    assert data["verdict"] == "DEAD_FRAGMENTED"
    assert data["non_vacuity_gate"]["pass"] is True
    rows = data["death_events"]
    assert len(rows) == 3
    assert [r["B"] for r in rows] == [98, 68, 56]
    assert [r["after_shape"]["live_shape"]["fragment_count_B"] for r in rows] == [2, 3, 4]
    assert [r["after_shape"]["live_shape"]["maxrun_B"] for r in rows] == [9, 4, 4]


def check_omega_extensions() -> None:
    p1e7 = load("PHASE0_OMEGA_LADDER_p1e7_2026-07-03.json")
    p1e8 = load("PHASE0_OMEGA_LADDER_p1e8_2026-07-03.json")
    assert p1e7["p"] == 10000019
    assert p1e7["max_omega_ge5"]["omega_ge5"] == 10
    assert p1e8["p"] == 100000007
    assert p1e8["max_omega_ge5"]["omega_ge5"] == 11


def main() -> None:
    check_mcolor_ladder()
    check_l3_deathrow()
    check_omega_extensions()
    print("OK: public Route-A obstruction packet JSON checks passed.")


if __name__ == "__main__":
    main()
