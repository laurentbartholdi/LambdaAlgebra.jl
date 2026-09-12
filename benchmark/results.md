# Measurements and remaining target

Measured on 2026-09-11 with Julia 1.11.7, Linux x86-64, one Julia thread,
`p=3`, `dimension=3`, and the default cache limit of 50,000 entries.
These are measurements of the implementation based on upstream `c8f3198`.

| Total-degree cutoff | Time to extend from previous row | Basis entries, all suffix contexts | Explicit boundary rows | Positive-degree homology classes |
| ---: | ---: | ---: | ---: | ---: |
| 40 | 0.094 s | 689 | 207 | 45 |
| 60 | 0.095 s | 2,927 | 940 | 117 |
| 80 | 0.827 s | 8,162 | 2,756 | 240 |
| 100 | 8.672 s | 19,558 | 6,699 | 466 |
| 120 | 184.798 s | 41,450 | 14,339 | 829 |

Constructor time, some compilation, and homology-copy time are excluded.
These basis counts are not byte/RAM measurements and do not include sparse
polynomials in the arithmetic caches. Timings are workload- and cache-dependent;
other experiments were also running during part of this session.

An earlier unpruned upstream run at total degree 64 retained 230,706 unstable
basis entries; the contextual prototype needed 3,668 entries **including its
auxiliary contexts**, and returned the identical homology signature.

## Validation

The optimized test suite passes 9,984 checks, including the explicit v0
tag-propagation counterexample. Independent full-matrix oracles
cover p=3 at sphere dimensions 1, 3, 5 (total cutoffs 24, 28, 24), and p=5 at
sphere dimensions 3 and 7 (cutoff 32). Comparisons with the usual lambda-mu
implementation cover S³ at p=3 through total degree 64, plus other odd spheres
and p=5 through 40. Prototype comparisons with the usual algebra also agreed
through 80; this larger comparison is not part of the routine test suite.

The large-run table above reports computed homology counts, not independent
full-matrix certification at total degree 120. The low-degree independent
oracles and the contextual cancellation argument validate the algorithm;
the ordinary algebra comparison is explicitly limited to the ranges listed
and to the additional benchmark comparison below.

## Periodic versus usual lambda-mu

Measured with Julia 1.11.7, one thread, p=3 and S³, after warming both
implementations at total degree 24. Each row starts a fresh algebra and
includes construction, Curtis reduction, truncation where applicable, and
`homology(A)`. All μ-degrees within the total cutoff are included. These
are single wall-clock measurements, not statistical estimates; they include
garbage collection and differ from the incremental timings above.

| Total-degree cutoff | Usual lambda-mu | Periodic | Positive-degree homology classes |
| ---: | ---: | ---: | ---: |
| 40 | 1.125 s | 1.834 s | 45 |
| 60 | 1.446 s | 2.303 s | 117 |
| 80 | 1.919 s | 3.321 s | 240 |
| 100 | 6.575 s | 12.508 s | 466 |

The positive-degree homology signatures agree in every tridegree at each
cutoff, extending the integrated comparison through total degree 100.
The usual implementation is faster in these runs. This is a comparison
of the current implementations, not an intrinsic complexity claim about
the two algebraic models. Explicit cycle reconstruction is not timed.

At cutoff 100, the periodic engine retains 8,806 `(μ,t,D)` contexts with
19,558 basis entries in total, plus its public basis and arithmetic caches.
These contexts are individual degree diagonals at required sphere bounds,
not 8,806 complete unstable algebras. The usual public basis has 999 entries
after truncation; this count excludes the stable entries discarded during
construction and is not a peak-memory comparison.

Reproduce with:

```sh
julia --project=. benchmark/compare_variants.jl 40 60 80 100
```

## Optimized implementation

The optimization keeps the public API and the contextual Curtis matching.
It adds interned words with compressed generator runs, direct binomial
formulas for products and differentials of v-powers, tighter saturation of
sphere bounds, two-generation disposable caches, and collection of temporary
word nodes. Construction reuses each Adem coefficient instead of recomputing
it for every generator pair. Expanded preimages and eager reduction of row
tails were tested and rejected because they were slower. An experimental
v₀-module computation was not adopted; it is not part of this implementation.

Fresh, sequential runs with the same Julia version, prime, sphere, warm-up
and complete total-degree cutoffs as above gave:

| Total-degree cutoff | Usual lambda-mu | Optimized periodic | Periodic allocated bytes | Positive-degree classes |
| ---: | ---: | ---: | ---: | ---: |
| 80 | 1.921 s | 0.514 s | 135,122,816 | 240 |
| 100 | 5.408 s | 3.408 s | 502,967,880 | 466 |
| 120 | 52.829 s | 63.585 s | 4,588,957,856 | 829 |

Each row includes construction and `homology(A)`, but not reconstruction of
cycle representatives. Allocated bytes are cumulative allocation traffic,
not peak RAM. The benchmark checks equality of the entire positive-degree
homology signature at every cutoff; the independent full-matrix oracle is
still limited to the smaller test ranges above.

At cutoff 100, periodic time fell from 12.508 s to 3.408 s (3.7 times faster)
and allocations from 3,963,982,704 to 502,967,880 bytes (7.9 times less).
The optimized periodic implementation is faster than the usual one at 80
and 100, but remains slower at 120. Intermediate boundary reconstruction
still grows rapidly; fewer basis monomials do not remove that cost.

The tighter sphere bound reduces auxiliary contexts at cutoff 100 from
8,806 to 2,916, and their basis entries from 19,558 to 16,572. In an incremental
run through 120, the collected word pool ended with 790,653 live nodes in
1,798,578 allocated slots; without collection and the binomial power formulas,
the compact-word prototype retained 7,309,513 nodes. These word counts are
not total-memory measurements.

The final suite passes 9,984 checks, including direct comparisons of the
compact differential with the public differential, p=3 and p=5 power formulas,
the user's unstable tag example, sphere-context saturation, and cache/word
collection. Run the comparison with:

```sh
julia --project=. benchmark/compare_variants.jl 80 100 120
```

## Dimension 500 is not yet achieved

A run requesting checkpoints through total degree 500 completed 120 and was
stopped during the next checkpoint, 140. No completed result at 140 or 500 is
claimed. Moreover, `top_degree=500` bounds total degree `μ+λ+top`; it would not
by itself establish completeness through topological degree 500.

The new redundancy criterion solves the unstable suffix problem and reduces
basis storage substantially. Completion of implicit boundaries still creates
large sparse expressions and grows rapidly beyond total degree 100. Caching
normal forms improves this but does not remove the growth. Reaching the stated
500 target needs a further improvement to that reduction stage, rather than
just increasing the generator array or omitting more monomials without proof.
The optimization above improves the achieved range's time and allocation
costs, but does not complete the dimension-500 target. Its larger incremental
run was stopped during cutoff 140 after completing 120.

Run the same benchmark with:

```sh
julia --project=. benchmark/unstable_periodic.jl 100
```

Passing `500` requests that total-degree target; this is a long-running
experiment, not a verified performance promise. `curtis_stats` and the flushed
checkpoint output make the achieved range visible.
