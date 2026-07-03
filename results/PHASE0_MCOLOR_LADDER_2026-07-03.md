# Phase-0 m_color Ladder Audit

Date: 2026-07-03
Status: COMPLETE for static Segment-A hit-table proxy.

This audit checks whether omega-maximal witness rows also have many active
GCD colors. It is not a row-DP replay: true DP `m_color(y)` depends on the
reachable reservoir `L_y(W)`. The count below is the static active-color
proxy over good B-lanes and `D*`; unreachable tokens could only reduce the
certified dynamic count.

Definitions used:

```text
D* = {-26,-22,-20,-10,-8,4,8,10,14}
x_c(B) = 86 p^2 + B p + 51
good B-lanes: B even and B == p (mod 3), starting at B >= 44, spaced by 6
ell active iff ell | x_c(B)+d for at least one audited (B,d)
```

## Results

| p | omega_max | m_color_static W=6 | m_color_static W=40 | m_color_static W=160 | tracks? | construction note |
|---:|---:|---:|---:|---:|---|---|
| 26267 | 7 | 6 | 7 | 7 | W160 yes | static B-window extension, not DP L_y(W) |
| 99961 | 7 | 6 | 7 | 7 | W160 yes | static B-window extension, not DP L_y(W) |
| 950039 | 8 | 4 | 8 | 8 | W160 yes | static B-window extension, not DP L_y(W) |
| 10000019 | 10 | 9 | 10 | 10 | W160 yes | static B-window extension, not DP L_y(W) |
| 100000007 | 11 | 9 | 11 | 11 | W160 yes | static B-window extension, not DP L_y(W) |

## Active Prime Detail

### p=26267, n=690022333

- omega_ge5: `7`
- factors_ge5: `7, 11, 13, 17, 23, 41, 43`
- W=6: m_color_static=`6`, active=`[7, 11, 13, 17, 23, 41]`, inactive=`[43]`, B=`44..74`
- W=40: m_color_static=`7`, active=`[7, 11, 13, 17, 23, 41, 43]`, inactive=`[]`, B=`44..278`
- W=160: m_color_static=`7`, active=`[7, 11, 13, 17, 23, 41, 43]`, inactive=`[]`, B=`44..998`

### p=99961, n=9992202220

- omega_ge5: `7`
- factors_ge5: `5, 7, 11, 13, 19, 109, 241`
- W=6: m_color_static=`6`, active=`[5, 7, 11, 13, 19, 241]`, inactive=`[109]`, B=`46..76`
- W=40: m_color_static=`7`, active=`[5, 7, 11, 13, 19, 109, 241]`, inactive=`[]`, B=`46..280`
- W=160: m_color_static=`7`, active=`[5, 7, 11, 13, 19, 109, 241]`, inactive=`[]`, B=`46..1000`

### p=950039, n=902574201802

- omega_ge5: `8`
- factors_ge5: `7, 13, 17, 19, 29, 37, 41, 349`
- W=6: m_color_static=`4`, active=`[7, 13, 17, 19]`, inactive=`[29, 37, 41, 349]`, B=`44..74`
- W=40: m_color_static=`8`, active=`[7, 13, 17, 19, 29, 37, 41, 349]`, inactive=`[]`, B=`44..278`
- W=160: m_color_static=`8`, active=`[7, 13, 17, 19, 29, 37, 41, 349]`, inactive=`[]`, B=`44..998`

### p=10000019, n=100000520554935

- omega_ge5: `10`
- factors_ge5: `5, 7, 11, 13, 17, 19, 37, 47, 71, 167`
- W=6: m_color_static=`9`, active=`[5, 7, 11, 13, 17, 19, 37, 71, 167]`, inactive=`[47]`, B=`44..74`
- W=40: m_color_static=`10`, active=`[5, 7, 11, 13, 17, 19, 37, 47, 71, 167]`, inactive=`[]`, B=`44..278`
- W=160: m_color_static=`10`, active=`[5, 7, 11, 13, 17, 19, 37, 47, 71, 167]`, inactive=`[]`, B=`44..998`

### p=100000007, n=10000002113691590

- omega_ge5: `11`
- factors_ge5: `5, 7, 11, 13, 19, 31, 37, 47, 59, 61, 271`
- W=6: m_color_static=`9`, active=`[5, 7, 11, 13, 19, 31, 37, 47, 59]`, inactive=`[61, 271]`, B=`44..74`
- W=40: m_color_static=`11`, active=`[5, 7, 11, 13, 19, 31, 37, 47, 59, 61, 271]`, inactive=`[]`, B=`44..278`
- W=160: m_color_static=`11`, active=`[5, 7, 11, 13, 19, 31, 37, 47, 59, 61, 271]`, inactive=`[]`, B=`44..998`

## Interpretation

- W=6: partial growth, slower than omega.
- W=40: tracks omega exactly on all witness rows.
- W=160: tracks omega exactly on all witness rows.

The parity claim is therefore not an automatic consequence of omega alone. On the old six-lane core, the large-p witnesses do not track omega. On the wider static windows, especially W=160, the active-color count does track omega on the measured witnesses. This supports a parity obstruction for full/static wide-window certificates, but it does not by itself refute the DP-reachable weakened L3' route.

## Source Artifacts

- `docs\PHASE0_Q1_EXACT_edge052_p26267_q26293_2026-07-02.json`
- `docs\PHASE0_Q1_EXACT_p99961_q99971_2026-07-02.json`
- `docs\PHASE0_Q1_EXACT_p950039_q950041_2026-07-02.json`
- `docs\PHASE0_OMEGA_LADDER_p1e7_2026-07-03.json`
- `docs\PHASE0_OMEGA_LADDER_p1e8_2026-07-03.json`
