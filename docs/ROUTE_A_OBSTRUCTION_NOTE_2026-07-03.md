# Route A Obstruction Note

Date: 2026-07-03
Status: BANKABLE NOTE DRAFT. This is not a proof of Erdos #1212. It records the
reduction, the obstruction to the static full-window elementary certificate,
the negative test of the index-interval L3' escape, the remaining general
DP-discrepancy wall, the shared quantifier barrier with Route S, and the one
credible non-sieve future direction.

## Abstract

The Segment-A row-DP construction for Erdos #1212 empirically maintains a deep
reachable reservoir, but the intended static full-window elementary certificate
for Scaled-Window Reservoir Survival is obstructed. Exact omega ladders up to
`p=100000007`, combined with the active-colour ladder audit, show that at
scaled static windows the active-colour count equals omega on the measured
witness rows: `m_color_static=omega=7,7,8,10,11` for W=40 and W=160. Thus
`2^m_color_static/(log p)^2` moves from the scale `1` to `6`, so the
Legendre/Jacobsthal full-window error term exceeds the target scaled window.
The sharper explanation has two parts: the static lower-bound certificate is
below the linear-sieve parity threshold `s=2`, and the interval length is
polylogarithmic in `N~p^2`, far below the pointwise short-interval scale reached
by known methods. The subsequent L3' death-row Q2 test also came back negative
for the tractable index-interval escape: on non-vacuous death rows the
DP-surviving lanes fragmented rather than forming one continuous thread. This
does not refute the existence of a path, but it pushes Route A back to the
general DP-discrepancy problem. Route S's density-1 strip theorem hits the same
all-large-column quantifier barrier. A percolation route based on random
coprime colouring is a genuine alternative framework, but current literature
gives almost-sure random-model cluster results, not deterministic
prime-prime-deleted #1212 survival.

## 1. Reduction Target

The Route-A construction attempts to build an infinite 4-adjacent path through
successive prime-square segments `[p^2,q^2]`, where `q=nextprime(p)`, using
centres

```text
x_c(B) = 86 p^2 + Bp + 51.
```

The proof obligation reduces to a pointwise all-large-p survival statement:

**Scaled-Window Reservoir Survival (target).** For every sufficiently large
prime-square segment and every row in that segment, the row-DP reachable
reservoir contains enough allowed offsets in a window of size
`C (log p)^2` to continue the chain and hand off to the next segment.

Empirically this reservoir is healthy on tested edges. The obstruction is not
the phenomenon; it is the static full-window elementary certificate for the
phenomenon.

## 2. Omega Ladder Evidence

| p | q | q-p | rows | (log p)^2 | log(p^2)/loglog(p^2) | omega_max_ge5 | 2^omega/(log p)^2 | witness_n_factors_ge5 |
|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 26267 | 26293 | 26 | 1366560 | 103.552374 | 6.754358 | 7 | 1.236089 | 7, 11, 13, 17, 23, 41, 43 |
| 99961 | 99971 | 10 | 1999320 | 132.538471 | 7.340812 | 7 | 0.965757 | 5, 7, 11, 13, 19, 109, 241 |
| 950039 | 950041 | 2 | 3800160 | 189.454807 | 8.303671 | 8 | 1.351246 | 7, 13, 17, 19, 29, 37, 41, 349 |
| 10000019 | 10000079 | 60 | 1200005880 | 259.793069 | 9.281705 | 10 | 3.941599 | 5, 7, 11, 13, 17, 19, 37, 47, 71, 167 |
| 100000007 | 100000037 | 30 | 6000001320 | 339.321482 | 10.214924 | 11 | 6.035574 | 5, 7, 11, 13, 19, 31, 37, 47, 59, 61, 271 |

Source summary: `docs/PHASE0_OMEGA_LADDER_SUMMARY_2026-07-03.md`.

The old static elementary route needs `2^m_color_static` to be negligible
compared with the scaled window. After the active-colour audit, the wider
static windows have `m_color_static=omega` on the measured ladder, so the omega
trend is the relevant trend for that certificate.

There is one nonmonotone dip at `p=99961`, where `(log p)^2` grows while
`omega_max_ge5` remains temporarily fixed at 7. This does not affect the
larger-scale trend from `1.236` to `6.036`.

## 2b. Active-Colour Ladder Audit

The omega ladder alone is not enough: earlier Segment-A work explicitly showed
that `m_color(y)` need not equal `omega(y+1)`. The relevant audit therefore
counted active GCD colours on the omega-maximal witness rows, using the
Segment-A static hit table

```text
D* = {-26,-22,-20,-10,-8,4,8,10,14}
x_c(B) = 86 p^2 + Bp + 51
good B-lanes: B even and B == p (mod 3)
ell active iff ell | x_c(B)+d for some audited (B,d)
```

The output is `docs/PHASE0_MCOLOR_LADDER_2026-07-03.md`.

| p | omega_max | m_color_static W=6 | m_color_static W=40 | m_color_static W=160 |
|---:|---:|---:|---:|---:|
| 26267 | 7 | 6 | 7 | 7 |
| 99961 | 7 | 6 | 7 | 7 |
| 950039 | 8 | 4 | 8 | 8 |
| 10000019 | 10 | 9 | 10 | 10 |
| 100000007 | 11 | 9 | 11 | 11 |

This resolves the first audit issue. The parity obstruction is not an automatic
consequence of omega alone: on the six-lane core, `m_color_static` is smaller
and nonmonotone. But at the wider static windows relevant to the scaled
reservoir test, W=40 and W=160 already activate every prime factor on the
measured ladder, so `m_color_static=omega` and `2^m_color_static=2^omega`.
For these measured witness rows and the anchored static good-B windows tested,
the static active-colour count is monotone in the B-window width and bounded by
omega. Hence once W=160 already gives `m_color_static=omega`, wider anchored
static B-windows in this setup, including any anchored window with at least
`(log p)^2 ~= 339` lanes at `p=100000007`, also have
`m_color_static=omega`.

This also exposes the narrow-window / wide-window tension. A narrow core can
keep `m_color` below omega, but it is too thin to serve as the scaled reservoir
supply. Widening the static window supplies enough tokens, but it also activates
the full omega bundle and triggers the parity obstruction. This is the precise
failure mode of the static full-window elementary certificate.

## 3. Obstruction for the Elementary Lower Sieve

The static full-window certificate asks for a lower-bound sieve in a short
window

```text
W ~= (log p)^2.
```

On omega-maximal rows, the largest active small prime `z` is large enough that

```text
s = log W / log z < 2.
```

For edge052's exact witness, `z=43`, so `z^2=1849` while `W~=103.55`; hence
`s~=1.23`. This is below the linear lower-sieve threshold. In the `s<2` region,
the lower-bound sieve function is trivial, so Rosser-Iwaniec / Fundamental
Lemma style lower bounds cannot certify survival. Selberg `Lambda^2` is an
upper-sieve tool and does not supply the missing lower bound.

Thus the obstruction is structural for the intended static full-window
certificate. Standard linear lower-sieve technology is in the `s<2`
zero-lower-bound regime for this object. Fixing constants in Legendre's
triangle-inequality error is not enough; one would need new pointwise Type-II or
bilinear information in a single polylogarithmic interval.

This statement is deliberately scoped. It does not rule out a non-sieve method,
an exact finite CRT method, or a DP argument that works directly with the
reachable set `L_y(W)` instead of the full static `B x D*` window.

## 3b. Second Obstruction: Polylogarithmic Window Scale

There is an independent scale obstruction. The row-DP window has length

```text
W = (log p)^2 ~= (1/4)(log N)^2, where N ~= p^2.
```

This is a polylogarithmic interval, not an interval of length `N^theta` for a
fixed positive `theta`. The pointwise short-interval tools available in the
literature are calibrated on power-length intervals. For primes, the classical
Baker-Harman-Pintz pointwise result gives primes in intervals of length
`x^0.525` for sufficiently large `x`, still a power scale and vastly longer
than `(log x)^2`. Even RH-scale pointwise prime-gap technology would remain
around `x^{1/2}` up to logarithmic factors, not polylogarithmic.

Matomaki-Radziwill style results are much closer in interval length, but the
quantifier is almost-all / average, not pointwise for every adversarial segment.
That is the same quantifier mismatch as Route S: average short-interval control
does not certify every block in a deterministic infinite path.

Therefore the intended static full-window elementary certificate faces two
separate obstructions: the standard lower sieve is in the parity dead zone, and
known pointwise short-interval technology does not reach the polylogarithmic
Cramer-scale window. Either obstruction is enough to block that old
full-window certificate.

## 3c. The Precise Frontier Gap

The frontier is best described on three axes:

```text
pointwise for every N; polylogarithmic length (log N)^2; Type-II / lower-bound input.
```

Current technology reaches these axes separately, but not simultaneously.

| Axis retained | Result | Length scale | Quantifier | Direction |
|---|---|---:|---|---|
| Pointwise prime asymptotic | Guth-Maynard, arXiv:2405.20552 | `N^(17/30+o(1))` | every `N` | lower/asymptotic |
| Pointwise prime existence | Baker-Harman-Pintz | `N^0.525` | every `N` | lower |
| Pointwise almost-prime existence | Campbell, arXiv:2603.10356 | `N^(1/2)` scale between consecutive squares | every `N` | lower for `Omega<=3` |
| Pointwise Type-II / nilsequence control | MSTT I, arXiv:2204.03754 | power scale, e.g. `N^(5/8)` regime | every interval | cancellation / upper-control input |
| Polylog almost-primes | Matomaki, arXiv:2012.11565; Matomaki-Teravainen, arXiv:2207.05038 | `(log N) h`, and almost-all `(log N)^2.1`-scale variants | almost all intervals | lower for `P_2` |
| Higher uniformity at short scale | MSTT II, arXiv:2411.05770 | short intervals | almost all intervals | averaged higher-uniformity |

The pointwise lower-bound wall is still at power scale, around `N^(1/2)` or
worse. The static full-window Route-A certificate asks for control at
`(log p)^2 ~= (1/4)(log N)^2`. Bridging from power-scale `N^(1/2)` to
Cramer-scale `(log N)^2` while keeping the pointwise for-every-`N` quantifier
is beyond current technology.

The MSTT I "all intervals" input is genuinely pointwise, but it has the wrong
shape for this proof: it operates at power-length scale and supplies
cancellation / upper-control information, not a lower bound certifying a
surviving reservoir. Conversely, the polylog-length almost-prime literature
obtains lower bounds only after averaging over `x`; that averaging is exactly
the loss of pointwise control that Route A and the actual infinite path cannot
afford.

Campbell's 2026 result is closest in geometry: for every `n`, the interval
`(n^2,(n+1)^2)` contains an integer with at most three prime factors, counted
with multiplicity. But this is still a length `~N^(1/2)` statement. Pushing
pointwise lower-bound sieve technology from consecutive-square length down to
polylogarithmic length would require a new input; getting actual primes at the
consecutive-square scale is already Legendre's conjecture.

Thus the claim is not "impossible forever". The precise claim is:

```text
No current technique bridges power-scale to Cramer-scale pointwise lower bounds
in the way Route A requires.
```

## 3d. Parity vs Discrepancy: Two Objects, Two Walls

An earlier internal synthesis attributed the true bottleneck to a large-sieve
discrepancy problem on the DP-reachable set and flagged an earlier broad
"parity barrier" phrasing as a misattribution. The audit reconciles this rather
than erasing it.

There are two distinct obstructions for two distinct objects:

```text
static full-window certificate:
  object = B x D*
  evidence = m_color_static = omega on W=40 and W=160 ladder witnesses
  obstruction = parity / lower-sieve dead zone s<2

DP-reachable weakened route L3':
  object = L_y(W) subset B x D*
  evidence = death-row Q2 is non-vacuous and index-interval survivors fragment
  obstruction class = index-interval escape fails; general L3 falls back to
                      finite-window CRT discrepancy / large-sieve control
```

Thus the earlier "parity misattributed" warning remains correct for the
DP-reachable object. Parity is correctly attributed to the static full-window
certificate after the `m_color` audit. No single wall explains all of Route A:
there are two walls for two different proof objects.

## 3e. L3' Death-Row Q2: Index-Interval Escape Fails

The open question after the static audit was whether the DP-reachable set could
avoid the parity obstruction by supplying an honest reachable sub-interval. This
was the L3' weakened route: instead of certifying the whole static `B x D*`
window, try to show that surviving DP lanes remain interval-like through actual
lane deaths.

The non-vacuous death-row Q2 test was run on the known death-bearing edge052
regime:

```text
p = 26267, q = 26293
W_lanes = 11
B = 44..104, good-class lanes spaced by 6
bank = DEFAULT_D
radius = 50, width = 30
```

Artifact:
`docs/L3PRIME_PHASE0_DEATHROW_Q2_edge052_W11_2026-07-03.md`.

The run was intentionally stopped after the third death (`rows_done=702591` of
`1366560`, `scan_complete_interval=false`), because that already met the
non-vacuity gate for this diagnostic:

```text
distinct interior death rows >= 3
observed separators >= 2
```

Results:

| dead B | row offset | live fragments after death | live MAXRUN | live GAPRATIO |
|---:|---:|---:|---:|---:|
| 98 | 14800 | 2 | 9 | 0.909091 |
| 68 | 186856 | 3 | 4 | 0.818182 |
| 56 | 702591 | 4 | 4 | 0.727273 |

Verdict: `DEAD_FRAGMENTED`.

This rules out the tested edge052/W=11 index-interval invariant as stated: the
surviving lane set does not remain one continuous thread through deaths. The
mechanism is also the expected one. Lane deaths are absorbing in this model;
when they occur at separated interior B-lanes, they create permanent gaps in
the index set, so the number of fragments can only grow.

This does not prove that the underlying visible path does not exist, and it
does not prove that the full general DP cannot route through a fragmented
survivor set. It does show that the tractable parity-avoiding sub-interval
certificate fails in this tested configuration. What remains is the original
general L3 problem: finite-window DP dynamics / CRT discrepancy without a cited
theorem.

## 4. Shared Barrier with Route S

Route S can plausibly produce a density-1 single-block strip connectivity
statement. It cannot imply Erdos #1212, because an infinite 4-adjacent path with
`x -> infinity` must intersect every sufficiently large vertical line `{x=N}`.
Therefore it needs an all-large-x pointwise crossing property, not an
almost-all-x exceptional-set theorem.

This is recorded separately in
`docs/PHASE0_QUANTIFIER_GEOMETRY_DETERMINATION_2026-07-03.md`.

The two proposed certificates therefore lose the needed quantifier at the same
logical place: pointwise control in short intervals. Route A needs it inside
row-DP reservoir survival; Route S needs it to upgrade density-1 blocks to an
actual infinite path.

## 5. Conditional Theorem Layer

The following conditional statement remains valid as a reduction, but the omega
ladder suggests its elementary hypothesis is not credible at the target
strength.

**Conditional Route-A theorem.** If there is a pointwise all-large-p reservoir
survival theorem strong enough to guarantee a living subreservoir of size
`>> (log p)^2` throughout every segment `[p^2,q^2]`, then the Segment-A row-DP
construction yields an infinite visible 4-adjacent path after deleting
prime-prime vertices.

The condition is not supplied by current elementary sieves. It should be treated
as the open core, not as a routine lemma.

## 6. Percolation Scoping

The one credible alternative framework found in the literature pass is random
coprime-colouring percolation:

- Martineau, "On coprime percolation, the visibility graphon, and the local
  limit of the GCD profile", arXiv:1804.06486.
- Le Fourn, Liu, Martineau, "Percolative properties of the random coprime
  colouring", arXiv:2509.08452.

The 2025 paper proves strong cluster results for the random coprime colouring;
in particular, for the standard `Z^d` graph with `d>=2`, it gives a unique
infinite white cluster and no infinite black cluster almost surely.

This does not solve #1212. The result is random-model / almost-sure, whereas
#1212 is deterministic and deletes the arithmetic set of prime-prime vertices.
A proof route would need a new robustness theorem showing that the visible
cluster survives this deletion in the actual lattice.

Scoping details: `docs/PHASE0_PERCOLATION_SCOPING_2026-07-03.md`.

## 6b. Percolation Route Relocates, Not Escapes

The percolation reformulation does not remove the all-large-column burden. In a
planar nearest-neighbour setting, a left-right open crossing of a rectangle is
blocked exactly by an 8-adjacent / king-move closed top-bottom dual crossing.
Thus a deterministic proof must still show, for every relevant large column
window, that the closed set

```text
non-visible cells union prime-prime cells union boundary/axis obstructions
```

does not contain such a top-bottom blocking path.

That statement is again pointwise in the horizontal coordinate. Soft density
input does not suffice: a set of very small density can still block by forming a
single vertical column. Sparsity of prime-prime cells also does not suffice,
because prime-prime cells are added to the already dense non-visible CRT
obstructions rather than acting alone.

The only cheap escape would be if the closed pattern were periodic modulo a
small fixed `Q`, reducing the all-large-column check to a finite period. A narrow
probe tested this door:

- artifact: `docs/PHASE0_PERCOLATION_PERIODICITY_PROBE_2026-07-03.md`
- windows: `512 x 512`, `y0=30000`
- positions: `x0 in {10000, 100000, 1000000, 10000000}`
- shifts: `Q in {30, 210}`
- result: `DOOR_CLOSED`

The closed pattern was not identical under either small primorial shift at any
tested position; changed-cell ratios were about `4.4%` to `7.3%`. The maximum
full closed/non-visible square varied from `2` to `3` over the tested positions.

This confirms the conceptual diagnosis: the percolation route relocates the
same pointwise short-interval arithmetic problem into a dual-blocking statement.
It does not escape it.

## References

- Tao, Terence. "Open question: The parity problem in sieve theory", 2007.
  https://terrytao.wordpress.com/2007/06/05/open-question-the-parity-problem-in-sieve-theory/
- "Parity problem (sieve theory)", Wikipedia, used only as a secondary pointer
  to the standard terminology and parity-sensitive sieve references.
  https://en.wikipedia.org/wiki/Parity_problem
- Matomaki, Kaisa; Radziwill, Maksym. "Multiplicative functions in short
  intervals II", arXiv:2007.04290.
  https://arxiv.org/abs/2007.04290
- Matomaki, Kaisa. "Almost primes in almost all very short intervals",
  arXiv:2012.11565.
  https://arxiv.org/abs/2012.11565
- Matomaki, Kaisa; Teravainen, Joni. "Almost primes in almost all short
  intervals II", arXiv:2207.05038.
  https://arxiv.org/abs/2207.05038
- Matomaki, Kaisa; Radziwill, Maksym; Shao, Xuancheng; Tao, Terence;
  Teravainen, Joni. "Higher uniformity of arithmetic functions in short
  intervals I", arXiv:2204.03754.
  https://arxiv.org/abs/2204.03754
- Matomaki, Kaisa; Radziwill, Maksym; Shao, Xuancheng; Tao, Terence;
  Teravainen, Joni. "Higher uniformity of arithmetic functions in short
  intervals II. Almost all intervals", arXiv:2411.05770.
  https://arxiv.org/abs/2411.05770
- Guth, Larry; Maynard, James. "New large value estimates for Dirichlet
  polynomials", arXiv:2405.20552.
  https://arxiv.org/abs/2405.20552
- Campbell, Peter. "On the Existence of Integers with at Most 3 Prime Factors
  Between Every Pair of Consecutive Squares", arXiv:2603.10356.
  https://arxiv.org/abs/2603.10356
- Baker, Harman, Pintz. "The difference between consecutive primes, II" (the
  `x^0.525` pointwise short-interval prime result). Bibliographic details to
  verify before venue submission; the exponent and statement are also cited in
  later short-interval-prime literature.
- Friedlander, Iwaniec. `Opera de Cribro`, for modern sieve-theory background
  including the beta sieve / fundamental lemma framework.
- Martineau, Sebastien. "On coprime percolation, the visibility graphon, and
  the local limit of the GCD profile", arXiv:1804.06486.
  https://arxiv.org/abs/1804.06486
- Le Fourn, Samuel; Liu, Mike; Martineau, Sebastien. "Percolative properties of
  the random coprime colouring", arXiv:2509.08452.
  https://arxiv.org/abs/2509.08452

## 7. Final Position

Bank the note with the narrowed framing. The static full-window elementary
certificate is shelved: it is parity-obstructed on the measured active-colour
ladder and also lies beyond the current pointwise Cramer-scale frontier. The
output is still useful: it prevents future work from spending months on a false
easy lemma and identifies the exact obstruction class for that certificate.

The DP-reachable weakened interval route L3' was then tested in a non-vacuous
death-bearing regime. Its index-interval version failed: the surviving lanes
fragmented on interior death rows. The honest status is:

```text
static full-window elementary route: parity-obstructed;
index-interval L3' route: tested negative / fragmented;
general DP L3 route: still possible in principle, but back to discrepancy wall.
```

Percolation was the only non-sieve alternative found in this pass. Sections 6
and 6b suggest that it relocates the same pointwise short-interval burden
through exact planar duality, rather than escaping it. The periodicity escape
hatch was empirically closed by the `DOOR_CLOSED` probe. We did not identify a
viable non-sieve route within this project.

The note therefore records negative evidence for the tractable Route-A
certificates explored here: static full-window sieve certificate,
index-interval L3', and the scoped percolation relocation. It should not be
advertised as a disproof of #1212 or of every conceivable DP route; the
remaining object is the open-ended general finite-window discrepancy problem.
