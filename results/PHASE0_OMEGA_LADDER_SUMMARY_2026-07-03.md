# Phase-0 Omega Ladder Summary

Date: 2026-07-03
Status: COMPLETE. This is the Etap 0 omega ladder summary for the Route-A shelve decision.

## Results

| p | q | q-p | rows | (log p)^2 | log(p^2)/loglog(p^2) | omega_max_ge5 | 2^omega | 2^omega/(log p)^2 | witness_n_factors_ge5 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 26267 | 26293 | 26 | 1366560 | 103.552374 | 6.754358 | 7 | 128 | 1.236089 | 7, 11, 13, 17, 23, 41, 43 |
| 99961 | 99971 | 10 | 1999320 | 132.538471 | 7.340812 | 7 | 128 | 0.965757 | 5, 7, 11, 13, 19, 109, 241 |
| 950039 | 950041 | 2 | 3800160 | 189.454807 | 8.303671 | 8 | 256 | 1.351246 | 7, 13, 17, 19, 29, 37, 41, 349 |
| 10000019 | 10000079 | 60 | 1200005880 | 259.793069 | 9.281705 | 10 | 1024 | 3.941599 | 5, 7, 11, 13, 17, 19, 37, 47, 71, 167 |
| 100000007 | 100000037 | 30 | 6000001320 | 339.321482 | 10.214924 | 11 | 2048 | 6.035574 | 5, 7, 11, 13, 19, 31, 37, 47, 59, 61, 271 |

## Trend

```text
p=26267     ratio=1.236089
p=99961     ratio=0.965757
p=950039    ratio=1.351246
p=10000019  ratio=3.941599
p=100000007 ratio=6.035574
```

The ratio `2^omega/(log p)^2` is not bounded near zero; it rises sharply by
`p~1e7` and `p~1e8`. This supports the shelve decision: the elementary
Legendre/Jacobsthal error term exceeds the scaled window, and the lower-sieve
route is parity-obstructed in the relevant pointwise short-interval regime.

The dip at `p=99961` is local: `(log p)^2` increases while `omega_max_ge5`
temporarily stays at 7. It does not change the larger trend from `1.236089` to
`6.035574`.

## Source Artifacts

- `docs/PHASE0_Q1_EXACT_edge052_p26267_q26293_2026-07-02.json`
- `docs/PHASE0_Q1_EXACT_p99961_q99971_2026-07-02.json`
- `docs/PHASE0_Q1_EXACT_p950039_q950041_2026-07-02.json`
- `docs/PHASE0_OMEGA_LADDER_p1e7_2026-07-03.json`
- `docs/PHASE0_OMEGA_LADDER_p1e8_2026-07-03.json`

## Watchdog Guard

The omega ladder runner `prime_square_phase0_omega_ladder.py` now has a
single-process watchdog thread: if `rows_done` does not increment within
`--no-progress-timeout-s` seconds, it atomically writes a progress JSON with
`status:"stalled"` and exits with code `2`. Verification run:

```text
python .\prime_square_phase0_omega_ladder.py --target 26267 --p 26267 --q 26293 --chunk-size 200000 --no-progress-timeout-s 2 --simulate-stall-after-chunks 1 --progress-json .\docs\_omega_watchdog_test_progress.json --output-json .\docs\_omega_watchdog_test.json --output-md .\docs\_omega_watchdog_test.md
```

Result: exit code `2`, progress file `docs/_omega_watchdog_test_progress.json`
with `status:"stalled"` and `rows_done:200000`.
