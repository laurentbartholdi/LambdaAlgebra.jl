# Unstable periodic Curtis reduction

This implementation is based on commit
`c8f319880c51eda2e140f6ce9504ce01cfc8e683`.

## Usage and scope

```julia
using LambdaAlgebra
A, lambda, v = periodic_algebra(p=3, dimension=3, top_degree=100)
H = homology(A)
LambdaAlgebra.curtis_stats(A)

# Increasing either cutoff fills any previously missing diagonals.
LambdaAlgebra.cache_basis!(A, 120)

# Unpruned reference computation, for small cutoffs only.
R, _, _ = periodic_algebra(p=3, dimension=3, top_degree=40, curtis=false)
```

The existing keyword `top_degree` actually bounds **total degree**
`μ + λ + top`, not `top` alone. A computation with `top_degree=500` is
therefore not a complete computation in topological degrees through 500.
The differential preserves total degree and μ. Each included diagonal is
computed completely, so homology in the included tridegrees is exact.
To interpret a sphere's homotopy degree, also account for the sphere shift
used by the package's printing routines.

This change makes finite, positive, odd sphere dimensions use contextual
Curtis tables. Stable computations retain the existing algorithm. Truncating
a periodic chain complex rebuilds its cancellations at the new sphere
bound. Truncate **before** taking homology; filtering a stable homology
basis cannot recover classes killed by sources that have become forbidden.

The generator capacity has been raised from 100 to 256 slots: the periodic
algebra has λ₁,…,λ₂₅₀, with up to six v-generators. Only v-generators whose
differentials fit this capacity are constructed. For p=3 this includes
v₀,…,v₅ and λ₁₂₅ of topological degree 499. Requests beyond the supported
periodic total-degree range fail explicitly.

## Why a single unstable suffix table is insufficient

Write an admissible monomial as

```math
m=v_{j_1}\cdots v_{j_a}\lambda_{i_1}\cdots\lambda_{i_b},
\qquad j_1\leq\cdots\leq j_a.
```

For `b>0`, the package's instability condition on an odd sphere of
dimension D is

```math
2i_1+1\leq D+2\sum_{r=1}^a p^{j_r}.
```

Pure v-monomials are also included. Removing a vᵢ from the left consequently
**increases** the sphere bound for the suffix by 2pⁱ. In particular, a suffix
forbidden on the original sphere can occur in an allowed full monomial.
Discarding it solely because it is forbidden on that sphere is incorrect.

For a pure lambda prefix λᵢ, admissibility gives a different suffix bound:
the next index is at most pi−1, so the suffix belongs to the pure lambda
complex with sphere bound 2pi−1. The prefix itself must satisfy 2i+1≤D.

We therefore use these two context transitions:

| Removed first letter | Suffix sphere bound | Further restriction |
| --- | --- | --- |
| vᵢ | D+2pⁱ | Its first v-index, if any, is at least i |
| λᵢ | 2pi−1 | The suffix contains only lambdas |

## Redundant monomials

A table is indexed by `(μ,t,D)`, where `t=λ+top`. Within it, entries are
ordered by the package's usual monomial lexicographic order. Suppose a
suffix entry x is paired with y by the Curtis algorithm **in the suffix
context just specified**. If both gx and gy are admissible, omit both
prefixed entries from the new table. An entry omitted recursively has the
same status: its pairing can be recovered by stripping the common prefix,
updating the context at each letter, until an explicitly stored tag is found.

Thus redundancy is a property of a monomial **and its suffix context**, not
of the suffix alone. Only the primitive tags and unpaired entries are stored
in the basis. A complete ordinary vector-space basis is never required.

### Justification

Use the usual Curtis leading-term extension property: if a completed tag
X has leading source x and normalized leading boundary y, then prepending
a generator g to an admissible source/target pair gives an implied tag
with leading source gx and leading boundary gy, with the sign changed for
a lambda prefix. This is the algebraic property underlying the existing
stable pruning algorithm; see Zhang, Theorem 3.2, cited below.

The additional issue here is ensuring that this cancellation takes place
inside the unstable complex. The contextual construction supplies exactly
that condition:

1. All monomials in the completed suffix X belong to its suffix context,
   because its elimination only uses sources in that complex.
2. For a vᵢ prefix, left multiplication lowers the instability bound by
   2pⁱ. Reordering the commuting v-prefix preserves its weight. The product
   vᵢX is therefore in the requested complex, including terms in X whose
   first v-index is smaller than i.
3. For a λᵢ prefix, the suffix is pure lambda and its first index is at
   most pi−1. Every term of λᵢX is admissible at the new join and has first
   lambda index i, which is allowed on the requested sphere.
4. The standard leading-term property gives an invertible pivot between
   gx and gy. Eliminating that source together with its boundary removes
   an acyclic pair. Lower terms are retained by reduction; they are not
   discarded as if the differential were a monomial map.

Induct on total degree. Every suffix context has smaller total degree,
including stripping v₀, which reduces μ by one. Within a diagonal the
algorithm proceeds by increasing lambda length, completing the incoming
map before using its target as a source. The inherited cancellations and
ordinary field elimination consequently preserve homology at every step.

For a fixed `t`, every first lambda index is at most `floor(t/(2p−2))`.
Bounds above `2floor(t/(2p−2))+1` describe the same complex. Canonicalizing
D at that value prevents an unbounded family of identical auxiliary tables.

The implementation also checks that every pivot produced during elimination
satisfies the current instability bound. A missing implied tag, a wrong
leading coefficient, or a failure of the expected triangular order raises
an error rather than silently dropping a term.

## Spectral-sequence interpretation

There are two compatible filtrations. Throughout, fix μ and t=λ+top,
so that the complex under consideration is finite dimensional. Write a
monomial as PL, with P a nondecreasing product of μ v-generators and L
an admissible pure lambda word. Set

```math
w(P)=\sum_j a_jp^j,\qquad
\tau(P)=\sum_j a_j(2p^j-2),\qquad P=\prod_j v_j^{a_j}.
```

First filter by the lexicographic order of the polynomial prefix P. The
Leibniz formula is

```math
d(PL)=d(P)L+P\,d(L).
```

Every term of d(P)L, after normalization, has polynomial prefix strictly
smaller than P: differentiation replaces some v_j by v_(j-1), and moving
the resulting lambda past later v-generators either preserves those later
indices or lowers them further. Hence the associated-graded differential
is just the pure-lambda differential. Its P-piece is

```math
P\,C^\lambda_{D+2w(P)}[t-\tau(P)],
```

where C^λ_B[u] denotes the pure-lambda complex with first index at most
(B−1)/2 and lambda internal degree u. In particular, the first page of
this coarser spectral sequence is

```math
E_1(P)\cong P\,H(C^\lambda_{D+2w(P)}[t-\tau(P)]).
```

The shifted bound depends on P; replacing all these complexes by C^λ_D
would already give the wrong associated graded. Higher differentials
incorporate the d(P)L term and corrections needed to lift cycles.

Next refine this filtration by the lambda word, using the full monomial
order employed by Curtis reduction. The usual lambda leading-term
triangularity, together with the strict decrease of polynomial prefixes
above, gives d(m) as a combination of smaller monomials. If all the allowed
monomials in the fixed diagonal are numbered m_1<...<m_N, then

```math
F_rC=\langle m_1,\ldots,m_r\rangle,\qquad dF_rC\subseteq F_{r-1}C.
```

Thus the refined associated-graded differential is zero, its E_1 page
has one symbol per allowed monomial, and its E_∞ page is the associated
graded of H(C_D). A Curtis tag x→y records completed representatives
X=x+(smaller terms) and dX=c y+(smaller terms), c≠0. It describes a
cancellation in this refined spectral sequence. With ordinal filtration
indices, its filtration drop is rank(x)−rank(y). The implementation jumps
directly between pivots; it does not enumerate these pages or store those
filtration drops in the package's `page` field. These auxiliary pages must
also not be confused with pages of the Adams spectral sequence.

The context transitions in the preceding section are how the same
filtered cancellation is reused after adding a prefix. It is essential
that the completed preimage be in the appropriate unstable complex,
not just that the target monomial be allowed.

### The first example on S³ at p=3

Here dλ₂=λ₁². The target λ₁² is allowed on S³ but λ₂ is not, so the
stable tag cannot kill that unstable class. After prefixing by v₀,
however, both v₀λ₂ and v₀λ₁² are allowed on S³ and

```math
d(v_0\lambda_2)=v_0\lambda_1^2.
```

Stripping v₀ sends the suffix calculation to S⁵, where λ₂→λ₁² is a
valid tag. The algorithm can therefore omit both prefixed monomials
without ever adding this pair to the S³ basis.

One sees the necessary correction to another source immediately:

```math
d(v_1\lambda_1)=v_0\lambda_1^2,
\qquad d(v_1\lambda_1-v_0\lambda_2)=0.
```

The omitted tag still participates in reconstructing this cycle. Omitting
a paired source/target means recovering their cancellation when needed,
not replacing both monomials by zero inside every expression.

### A tag that cannot be propagated from S³

At p=3, put

```math
x=v_1^3\lambda_2\lambda_1,\qquad
z=v_0^3\lambda_3\lambda_2\lambda_1,\qquad
a=v_0^3\lambda_5\lambda_1.
```

On S³ the Curtis table pairs x with z. More precisely, the raw differential
is `d(x)=z+v₁³λ₁³`; eliminating the larger term gives the completed source
`X=v₁³(λ₂λ₁+λ₁λ₂)`, with

```math
dX=v_0^3\lambda_3(\lambda_2\lambda_1+\lambda_1\lambda_2).
```

Both v₀x and v₀z are allowed on S³. Nevertheless, the S³ tag x→z must
not be used to omit v₀x: the smaller source v₀a becomes allowed and owns
the pivot v₀z. Indeed, a first appears on S⁵, whereas v₀a first appears
on S³. In the inequality `k≤n+w(P)`, a requires `5≤1+3`, which fails,
but v₀a requires `5≤1+4`, which holds.

The implementation removes v₀ and consults the **S⁵** suffix table.
In that table a owns z with pivot coefficient −1, and x survives.
Consequently v₀x is stored as an unpaired entry on S³. The omitted pivot
v₀z is recovered from v₀a, not v₀x. With `curtis=false`, these same
source and target entries are explicit and give the same answer.

The cycle returned by `cycle` for the leading monomial v₀x is

```math
C=v_0v_1^3(\lambda_2\lambda_1+\lambda_1\lambda_2)
 +v_0^4(\lambda_5\lambda_1+\lambda_4\lambda_2-\lambda_2\lambda_4)
 -v_0^3v_1\lambda_2\lambda_3.
```

Every term is allowed on S³, `dC=0`, and C represents the retained class
in `(μ,λ,top)=(4,2,22)`. Thus v₀a is a correction to the cycle led by
v₀x; it is not the differential target of v₀x (both have lambda length
two). The regression test checks the tags in both sphere contexts, this
exact cycle, and agreement with pruning disabled.

This example illustrates why a converse tag-propagation statement between
tables at the same fixed sphere bound fails. Prefixing by v₀ changes which
smaller sources can compete. The appropriate suffix table must already
include those sources before its cancellations can be reused.

## Arithmetic and storage

The private engine uses sparse dictionaries over GF(p), memoizes generator
multiplication in admissible form, and builds differentials recursively by
the Leibniz rule. It retains reduced rows for explicit tags. Implied rows
are reconstructed from suffix tags and cached separately. Prefixing a cycle
has zero Leibniz correction, so those rows can be lifted directly.

`cache_limit` bounds the **number of entries in each disposable cache**:
normalized products, suffix differentials, and implied rows. When a cache
fills, it is cleared and entries are recomputed as needed. Setting the limit
to zero disables those caches. This is not a bound on total RAM: explicit
tables/rows are retained, and one sparse row can have many terms.

`curtis_stats(A).stored` includes the auxiliary contexts. Public bases are
copies of their context tables, so `homology`, copying, and truncation cannot
delete the tags required to interpret omitted monomials. Cycle representatives
remain available through the existing `cycle(A, monomial)` interface.

## Comparison convention and validation

The current usual lambda-mu code excludes every nonempty word ending in μ.
It therefore omits the degree-zero tower. The periodic algebra includes
v₀ᵏ for every k; these are nonzero classes in `(μ,λ,top)=(k,0,0)`.
This difference is present in the upstream code and has not been hidden by
changing either complex. Comparisons with the usual code are made in
**positive topological degree**; the periodic tower is checked separately.

The test suite uses an independent enumeration by polynomial prefixes and
left-to-right lambda words, followed by full matrices and ordinary modular
Gaussian elimination. This oracle uses neither suffix tables nor pruning.
It also checks d²=0, compares the memoized differential with the public one,
checks completed cycles and their instability, compares with the usual
lambda-mu algebra at p=3 and p=5, and exercises incremental cutoffs, truncation,
copying, and disabled caches.

Performance measurements and the outstanding dimension-500 goal are recorded
in `benchmark/results.md`. There is no claim here that dimension 500 has been
computed or that a cutoff on μ gives a complete answer in every filtration.

## References

- Brayton Gray, *The periodic lambda algebra*, Fields Institute
  Communications 19 (1998), 93–101.
  [Paper](https://www.sas.rochester.edu/mth/sites/doug-ravenel/otherpapers/gray3.pdf).
- Jessie Zhang, *Periodicity in the periodic lambda algebra* (2012),
  §3, especially Theorem 3.2 for the stable leading-term extension property.
  [Paper](https://www.sas.rochester.edu/mth/sites/doug-ravenel/otherpapers/2012Zhang.pdf).
