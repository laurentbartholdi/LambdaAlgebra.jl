const LAMBDA = (:λ,)

index(g::Gen{2,LAMBDA}) = g.x-1
kind(g::Gen{2,LAMBDA}) = 'λ'

is_admissible_pair(g::Gen{2,LAMBDA},h::Gen{2,LAMBDA}) = 2*index(g)≥index(h)
is_admissible_product(g::Gen{2,LAMBDA},m::Monomial{2,LAMBDA}) = isempty(m) || is_admissible_pair(g,m[1])

# the first sphere dimension at which this appears
dimension(Λ::Algebra{2,LAMBDA},m::Monomial{2,LAMBDA}) = isempty(m) ? 0 : m[1].x

"""lambda2_algebra(p=2)

Construct the lambda algebra in characteristic 2
"""
function lambda2_algebra(;top_degree=-1, dimension=STABLE_DIMENSION)
    K = GF{2}

    A(k,j) = K(binomial(big(k-j-1),j))

    GEN = Gen{2,LAMBDA}
    NLAMBDA = NGEN-1
    
    λGen(i) = GEN(i+1)
    
    rules = Matrix{Union{Nothing,Vector{Tuple{K,GEN,GEN}}}}(undef,NGEN,NGEN)
    fill!(rules,nothing)
    for i=0:NLAMBDA, k=0:NLAMBDA-1-2i # λᵢλ₂ᵢ₊₁₊ₖ
        rule = []
        for j=0:k-1
            iszero(A(k,j)) || push!(rule,(A(k,j),λGen(i+k-j),λGen(2i+1+j)))
        end
        rules[λGen(i),λGen(2i+1+k)] = rule
    end
    
    diff = Vector{Vector{Tuple{K,GEN,GEN}}}(undef,NGEN)
    for k=1:NLAMBDA+1
        diff[λGen(k-1)] = [(A(k,j),λGen(k-j-1),λGen(j-1)) for j=1:k-1 if !iszero(A(k,j))]
    end
    
    degrees = Vector{Degree}(undef,NGEN)
    for i=0:NLAMBDA
        degrees[λGen(i)] = (μ=0,λ=1,top=i)
    end
                                            
    Λ = Algebra{2,LAMBDA}(rules,diff,degrees,Dict(),Ref(-1),Ref(1),Ref(dimension))
    # sanity check
    for g::GEN=1:NGEN, h::GEN=1:NGEN
        @assert is_admissible_pair(g,h) == (rules[g,h]==nothing)
    end

    cache_basis!(Λ,top_degree)
    if dimension≠STABLE_DIMENSION
        truncate!(Λ,dimension)
    end

    Λ, [AlgebraElem(Λ,Dict(Monomial(λGen(i))=>one(K))) for i=0:NLAMBDA]
end
