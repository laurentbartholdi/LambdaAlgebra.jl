const LAMBDAV = (:v,:μ)
const NV = 5

# 1:NGEN-NV-1 are λ, NGEN-NV:NGEN are v
index(g::Gen{p,LAMBDAV}) where p = g.x≤NGEN-NV-1 ? g.x : g.x-(NGEN-NV)
kind(g::Gen{p,LAMBDAV}) where p = g.x≤NGEN-NV-1 ? 'λ' : 'v'

is_admissible_pair(g::Gen{p,LAMBDAV},h::Gen{p,LAMBDAV}) where {p} = is_lambda(g) ? is_lambda(h) && index(h)≤p*index(g)-1 : is_lambda(h) || index(h)≥index(g)
is_admissible_product(g::Gen{p,LAMBDAV},m::Monomial{p,LAMBDAV}) where {p} = isempty(m) || is_admissible_pair(g,m[1])

# the first sphere dimension at which this appears
function dimension(Λ::Algebra{p,LAMBDAV},m::Monomial{p,LAMBDAV}) where {p}
    dim = 0
    for g=m
        if is_lambda(g)
            return max(0,2index(g)+1-dim)
        else
            dim += Λ.degree[g].top+2
        end
    end
    return 0
end

"""periodic_algebra(; p=3, top_degree=-1, μ_degree=-1,
                       dimension=STABLE_DIMENSION, curtis=true, cache_limit=50000)

Construct the periodic lambda algebra in odd characteristic `p`.
`top_degree` is the existing total-degree cutoff `μ+λ+top`.
Positive odd sphere dimensions use context-sensitive Curtis reduction.
For small reference computations, `curtis=false` disables pruning.
`cache_limit` bounds entries in each disposable arithmetic cache, not total RAM.
See `docs/unstable-periodic.md` for the criterion and comparison conventions.
"""
function periodic_algebra(;p=3, top_degree=-1, μ_degree=-1, dimension=STABLE_DIMENSION, curtis=true, cache_limit=50000)
    @assert isodd(p)
    K = GF{p}

    NLAMBDA = NGEN-NV-1
    max_v = min(NV, ndigits(NLAMBDA; base=p))
    ngen = NGEN-NV+max_v
    if dimension != STABLE_DIMENSION
        isodd(dimension) && dimension > 0 || throw(ArgumentError("the unstable periodic algebra requires a positive odd sphere dimension"))
    end
    
    # The same Adem coefficient is used for every possible first index.
    # Compute it once, rather than repeatedly evaluating large binomials.
    coefficients=zeros(K,NLAMBDA,NLAMBDA)
    for k=1:NLAMBDA, j=0:k-1
        coefficients[j+1,k]=K((-1)^(j+1)*binomial(big(p-1)*(k-j)-1,j))
    end
    A(k,j) = coefficients[j+1,k]

    GEN = Gen{p,LAMBDAV}
    λGen(i) = GEN(i)
    vGen(i) = GEN(i+NGEN-NV)
    
    rules = Matrix{Union{Nothing,Vector{Tuple{K,GEN,GEN}}}}(undef,ngen,ngen)
    fill!(rules,nothing)
    for i=0:max_v, j=i+1:max_v # all v's increasing
        rules[vGen(j),vGen(i)] = [(one(K),vGen(i),vGen(j))]
    end
    for i=1:NLAMBDA, j=0:max_v # move v's before λ's
        if j==0 || i+p^(j-1) > NLAMBDA
            rules[λGen(i),vGen(j)] = [(one(K),vGen(j),λGen(i))]
        else
            rules[λGen(i),vGen(j)] = [(one(K),vGen(j),λGen(i)),(one(K),vGen(j-1),λGen(i+p^(j-1)))]
        end
    end
    for i=1:NLAMBDA, k=0:NLAMBDA-p*i # λ's don't increase too much
        rule = Tuple{K,GEN,GEN}[]
        for j=0:k-1
            iszero(A(k,j)) || push!(rule,(A(k,j),λGen(i+k-j),λGen(p*i+j)))
        end
        rules[λGen(i),λGen(p*i+k)] = rule
    end
    
    diff = Vector{Vector{Tuple{K,GEN,GEN}}}(undef,ngen)
    for k=1:NLAMBDA
        diff[λGen(k)] = [(A(k,j),λGen(k-j),λGen(j)) for j=1:k-1 if !iszero(A(k,j))]
    end
    diff[vGen(0)] = [] # d(v₀) = 0
    for j=1:max_v
        diff[vGen(j)] = [(one(K),vGen(j-1),λGen(p^(j-1)))]
    end
    
    degrees = Vector{Degree}(undef,ngen)
    for i=1:NLAMBDA
        degrees[λGen(i)] = (μ=0,λ=1,top=(2p-2)*i-1)
    end
    for j=0:max_v
        degrees[vGen(j)] = (μ=1,λ=0,top=2*p^j-2)
    end
                                            
    Λ = Algebra{p,LAMBDAV}(rules,diff,degrees,dimension)
    # sanity check
    for g::GEN=1:ngen, h::GEN=1:ngen
        @assert is_admissible_pair(g,h) == (rules[g,h]==nothing)
    end

    if dimension != STABLE_DIMENSION
        Λ.curtis = PeriodicCurtis(Λ; prune=curtis, cache_limit)
    end
    cache_basis!(Λ,top_degree,μ_degree)
    
    Λ, [AlgebraElem(Λ,Monomial(λGen(i))) for i=1:NLAMBDA],[AlgebraElem(Λ,Monomial(vGen(i))) for i=0:max_v]
end
