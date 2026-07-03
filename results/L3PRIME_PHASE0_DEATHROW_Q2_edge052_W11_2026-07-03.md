# L3' Phase-0 Death-Row Q2

- generated_at: `2026-07-03T09:05:59`
- status: `complete`
- p: `26267`
- q: `26293`
- W_lanes: `11`
- B_range: `[44, 104]`
- rows_done: `702591` / `1366560`
- death_event_count: `3`
- distinct_death_rows: `3`
- interior_death_count: `3`
- max_after_live_fragments: `4`
- max_after_dstar_fragments: `4`
- non_vacuity_gate: `True`
- verdict: `DEAD_FRAGMENTED`

## Non-Vacuity Gate

```json
{
  "min_distinct_death_rows": 3,
  "min_interior_deaths": 3,
  "min_separators": 2,
  "observed_distinct_death_rows": 3,
  "observed_interior_deaths": 3,
  "observed_separators_from_live_fragments": 3,
  "pass": true
}
```

## Death Events

| B | death_y | row_offset | interior | after live count | after live MAXRUN | after live GAPRATIO | after live fragments | after D* MAXRUN | after D* fragments |
|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|
| 98 | 689970089 | 14800 | True | 10 | 9 | 0.909091 | 2 | 9 | 2 |
| 68 | 690142145 | 186856 | True | 9 | 4 | 0.818182 | 3 | 4 | 3 |
| 56 | 690657880 | 702591 | True | 8 | 4 | 0.727273 | 4 | 4 | 4 |

## Shape Extrema

```json
{
  "max_dstar_fragments": {
    "count": 8,
    "fragment_count_B": 4,
    "fragment_count_index": 4,
    "gapratio_B": 0.7272727272727273,
    "gapratio_index": 0.7272727272727273,
    "maxrun_B": 4,
    "maxrun_index": 4,
    "row_offset": 188597,
    "y": 690143886
  },
  "max_live_fragments": {
    "count": 8,
    "fragment_count_B": 4,
    "fragment_count_index": 4,
    "gapratio_B": 0.7272727272727273,
    "gapratio_index": 0.7272727272727273,
    "maxrun_B": 4,
    "maxrun_index": 4,
    "row_offset": 702591,
    "y": 690657880
  },
  "min_dstar_shape": {
    "count": 8,
    "fragment_count_B": 4,
    "fragment_count_index": 4,
    "gapratio_B": 0.7272727272727273,
    "gapratio_index": 0.7272727272727273,
    "maxrun_B": 4,
    "maxrun_index": 4,
    "row_offset": 188597,
    "y": 690143886
  },
  "min_live_shape": {
    "count": 8,
    "fragment_count_B": 4,
    "fragment_count_index": 4,
    "gapratio_B": 0.7272727272727273,
    "gapratio_index": 0.7272727272727273,
    "maxrun_B": 4,
    "maxrun_index": 4,
    "row_offset": 702591,
    "y": 690657880
  },
  "min_living_Dstar_pairs": 19,
  "min_living_Dstar_pairs_record": {
    "Dstar_counts": {
      "104": 2,
      "44": 1,
      "50": 0,
      "56": 1,
      "62": 2,
      "74": 3,
      "80": 2,
      "86": 3,
      "92": 5
    },
    "bank_sizes": {
      "104": 4,
      "44": 4,
      "50": 4,
      "56": 2,
      "62": 3,
      "74": 5,
      "80": 5,
      "86": 5,
      "92": 8
    },
    "dead_B": [
      68,
      98
    ],
    "dstar_nonempty_B": [
      44,
      56,
      62,
      74,
      80,
      86,
      92,
      104
    ],
    "dstar_shape": {
      "count": 8,
      "fragment_count_B": 4,
      "fragment_count_index": 4,
      "gapratio_B": 0.7272727272727273,
      "gapratio_index": 0.7272727272727273,
      "maxrun_B": 4,
      "maxrun_index": 4
    },
    "live_B": [
      44,
      50,
      56,
      62,
      74,
      80,
      86,
      92,
      104
    ],
    "live_shape": {
      "count": 9,
      "fragment_count_B": 3,
      "fragment_count_index": 3,
      "gapratio_B": 0.8181818181818182,
      "gapratio_index": 0.8181818181818182,
      "maxrun_B": 4,
      "maxrun_index": 4
    },
    "living_Dstar_nonempty_B_count": 8,
    "living_Dstar_pairs": 19,
    "row_offset": 208540,
    "y": 690163829
  }
}
```

## Interpretation

Death rows fragmented the surviving lane set. This rules out the tested edge052/W=11 index-interval invariant as stated; it does not rule out a fragmented general DP argument.
