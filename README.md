# Erdos #1212: Route-A Obstruction Notes

This repository contains a cleaned public packet from an exploratory,
AI-assisted attempt on Erdos problem #1212.

Problem #1212 asks for an infinite 4-adjacent path through visible lattice
points after deleting prime-prime vertices. We explored one constructive route:
a row-DP through successive prime-square blocks `[p^2, q^2]`,
`q = nextprime(p)`.

The outcome is **not** a proof or disproof of #1212. It is an obstruction note
for one family of constructive certificates.

## Main Document

- [`docs/ROUTE_A_OBSTRUCTION_NOTE_2026-07-03.md`](docs/ROUTE_A_OBSTRUCTION_NOTE_2026-07-03.md)

## What This Packet Claims

Allowed claims:

- We measured an omega ladder and an active-colour ladder on selected
  witness rows up to `p = 100000007`.
- In the measured wide static windows, `m_color_static` tracks the omega-maximal
  witnesses: `7, 7, 8, 10, 11`.
- The tested static full-window elementary certificate faces a parity/short
  interval obstruction.
- A non-vacuous L3' death-row diagnostic on edge052/W=11 showed fragmentation
  of the surviving lane set.
- These results obstruct the tested certificates, not the existence of a path.

Blocked claims:

- This is not a proof of Erdos #1212.
- This is not a disproof of Erdos #1212.
- This does not prove that every Route-A/DP approach is impossible.
- This does not claim novelty beyond the cited finite artifacts and discussion.
- This is not peer reviewed.

## Repository Layout

- `docs/` - cleaned notes and draft forum post.
- `results/` - selected JSON/MD finite artifacts used by the note.
- `src/` - small scripts for omega and active-colour audits.
- `experimental/` - heavier runner code for the L3' diagnostic and local
  dependencies. These are included for transparency, not as a polished package.
- `verify_results.py` - lightweight consistency checks for the published JSONs.

## Quick Verification

Install dependencies:

```bash
python -m pip install -r requirements.txt
```

Run the lightweight checks:

```bash
python verify_results.py
```

This verifies the expected finite summaries from the stored JSON files. It does
not rerun the expensive ladder scans.

## Reproducing Heavier Runs

The original exploratory workspace contained many intermediate scripts and
large worker outputs. This public repository intentionally excludes that
chaotic material.

The included scripts can be used as starting points for reproduction, but the
long omega ladder scans and L3' death-row runner are not a one-command polished
pipeline. Treat them as research artifacts.
