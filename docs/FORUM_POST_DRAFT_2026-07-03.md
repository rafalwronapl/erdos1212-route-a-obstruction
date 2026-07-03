# Forum Post Draft - Erdos #1212 Constructive Obstruction

Draft date: 2026-07-03
Status: draft for user review; do not publish without manual edit/link check.

## Draft

Obstructions we hit on one constructive approach to #1212 - feedback welcome

We explored a constructive route to Erdos #1212: an infinite 4-adjacent path
through visible lattice points after deleting prime-prime vertices. The route
uses a row-DP through successive prime-square blocks `[p^2,q^2]`,
`q=nextprime(p)`, maintaining a reservoir of admissible offsets in a window of
scale `(log p)^2`.

Empirically, the DP reservoir looks healthy on the tested finite blocks. But
the tractable certificates we hoped to use for the asymptotic proof appear to
hit real obstructions.

The issue is not just that rows with many prime factors exist. We audited the
active colours directly. On the omega-maximal witness rows in our ladder up to
`p=100000007`, the static active-colour count equals omega once the B-window is
wide enough for the scaled reservoir test:

```text
p=26267      omega=7   m_color_static(W=40)=7   m_color_static(W=160)=7
p=99961      omega=7   m_color_static(W=40)=7   m_color_static(W=160)=7
p=950039     omega=8   m_color_static(W=40)=8   m_color_static(W=160)=8
p=10000019   omega=10  m_color_static(W=40)=10  m_color_static(W=160)=10
p=100000007  omega=11  m_color_static(W=40)=11  m_color_static(W=160)=11
```

The old six-lane core does not behave this way (`m_color_static` is smaller and
nonmonotone there), so this is not the naive false shortcut
`m_color=omega`. The tension is that a narrow core can keep `m_color` low, but
is too thin to supply the scaled reservoir; widening the static window supplies
tokens but activates the full omega bundle. (`m_color_static` is nondecreasing
in the window width and bounded by omega, so equality already at `W=160`
carries over to the full scaled window `(log p)^2`, which reaches about `339`
at `p=10^8`.)

For the static full-window certificate, this puts the lower-bound sieve below
the parity threshold. At the edge052 witness, for example, the largest active
prime is `z=43` while `W~(log p)^2~103.55`, so

```text
s = log W / log z ~= 1.23 < 2.
```

In that `s<2` regime the linear lower sieve is trivial. Independently, the
window is polylogarithmic in `N~p^2`, i.e. Cramer-scale, whereas pointwise
short-interval lower-bound technology remains at power scale. We did not find
a result bridging power-scale to `(log N)^2` while keeping the pointwise
for-every-`N` quantifier.

We then tested the natural weakened escape (`L3'`): perhaps the DP-reachable
set, unlike the full static window, contains one honest surviving sub-interval
through death events. A non-vacuous death-row Q2 test on edge052/W=11 came back
negative. The run was intentionally stopped after the third death
(`rows_done=702591` of `1366560`), because it had already met the diagnostic
gate: 3 distinct interior death rows and 3 observed separators. After
successive absorbing lane deaths, the surviving lane set fragmented:

```text
B=98  row_offset=14800   live fragments=2  MAXRUN=9  GAPRATIO=0.909
B=68  row_offset=186856  live fragments=3  MAXRUN=4  GAPRATIO=0.818
B=56  row_offset=702591  live fragments=4  MAXRUN=4  GAPRATIO=0.727
```

This does not prove that no DP path exists. It does rule out the tested easy
index-interval certificate: the surviving reachable lanes do not form a single
continuous thread through deaths. What remains is the harder general
DP-discrepancy problem, not the clean sub-interval certificate we were hoping
for.

So the honest status is:

```text
static full-window elementary certificate: parity-obstructed;
index-interval L3' certificate: tested negative / fragmented;
general DP route: not refuted, but back to finite-window discrepancy.
```

We also scoped a percolation reformulation via coprime-colouring/planar duality.
It looks like that route relocates the same pointwise short-interval burden into
a dual blocking statement rather than escaping it, although this is only a
scoping result.

Caveat: this is exploratory, AI-assisted analysis and may contain errors. In
particular the `s<2` sieve setup has not been checked by a human expert, and the
"reduces to a discrepancy problem" framing is our own reading, not a cited
theorem. Corrections are very welcome.

We would genuinely welcome feedback on any of these points:

1. Is there a mistake in the `s<2` setup for the static full-window certificate?
2. Is there pointwise short-interval / Type-II technology that reaches a
   polylogarithmic lower-bound regime of this kind?
3. Is there a known way to control the fragmented DP-reachable set without
   reducing to a finite-window CRT discrepancy / large-sieve problem?

Full write-up:
https://github.com/rafalwronapl/erdos1212-route-a-obstruction/blob/main/docs/ROUTE_A_OBSTRUCTION_NOTE_2026-07-03.md
