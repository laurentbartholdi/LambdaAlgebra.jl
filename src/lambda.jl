const LAMBDA = (:λ,)

index(g::Gen{p,LAMBDA}) where p = g.x-(p==2)
kind(g::Gen{p,LAMBDA}) where p = 'λ'

#is_admissible_pair(g::Gen{2,LAMBDA},h::Gen{2,LAMBDA}) = 2*index(g)≥index(h)
is_admissible_pair(g::Gen{p,LAMBDA},h::Gen{p,LAMBDA}) where p = p*index(g)≥index(h)+(p≠2)

is_admissible_product(g::Gen{p,LAMBDA},m::Monomial{p,LAMBDA}) where p = isempty(m) || is_admissible_pair(g,m[1])

# the first sphere dimension at which this appears
dimension(Λ::Algebra{p,LAMBDA},m::Monomial{p,LAMBDA}) where p = isempty(m) ? 0 : 2m[1].x

"""lambda2_algebra(p=2)

Construct the lambda algebra in characteristic 2
"""
function lambda2_algebra(;top_degree=-1, μ_degree=-1, dimension=STABLE_DIMENSION)
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
                                            
    Λ = Algebra{2,LAMBDA}(rules,diff,degrees,dimension)
    # sanity check
    for g::GEN=1:NGEN, h::GEN=1:NGEN
        @assert is_admissible_pair(g,h) == (rules[g,h]==nothing)
    end

    cache_basis!(Λ,top_degree,μ_degree)
    if dimension≠STABLE_DIMENSION
        truncate!(Λ,dimension)
    end

    Λ, [AlgebraElem(Λ,Monomial(λGen(i))) for i=0:NLAMBDA]
end

"""pure_lambda_algebra(p=3)

Construct the pure lambda algebra in odd characteristic `p`
"""
function pure_lambda_algebra(;p=3, top_degree=-1, μ_degree=-1, dimension=STABLE_DIMENSION)
    p==2 && return lambda2_algebra(;top_degree,dimension)
    
    @assert isodd(p) # actually test isprime
    K = GF{p}

    A(k,j) = K((-1)^(j+1) * binomial(big(p-1)*(k-j)-1,j))
    B(k,j) = K((-1)^j * binomial(big(p-1)*(k-j),j))

    GEN = Gen{p,LAMBDA}
    NLAMBDA = NGEN
    
    λGen(i) = GEN(i)
    
    rules = Matrix{Union{Nothing,Vector{Tuple{K,GEN,GEN}}}}(undef,NGEN,NGEN)
    fill!(rules,nothing)
    for i=1:NLAMBDA, k=0:NLAMBDA-p*i # λᵢλⱼ
        rule = []
        for j=0:k-1
            iszero(A(k,j)) || push!(rule,(A(k,j),λGen(i+k-j),λGen(p*i+j)))
        end
        rules[λGen(i),λGen(p*i+k)] = rule
    end
    
    diff = Vector{Vector{Tuple{K,GEN,GEN}}}(undef,NGEN)
    for k=1:NLAMBDA
        diff[λGen(k)] = [(A(k,j),λGen(k-j),λGen(j)) for j=1:k-1 if !iszero(A(k,j))]
    end

    degrees = Vector{Degree}(undef,NGEN)
    for i=1:NLAMBDA
        degrees[λGen(i)] = (μ=0,λ=1,top=(2p-2)*i-1)
    end
                                            
    Λ = Algebra{p,LAMBDA}(rules,diff,degrees,dimension)
    # sanity check
    for g::GEN=1:NGEN, h::GEN=1:NGEN
        @assert is_admissible_pair(g,h) == (rules[g,h]==nothing)
    end

    cache_basis!(Λ,top_degree,μ_degree)
    if dimension≠STABLE_DIMENSION
#        truncate!(Λ,abs(dimension))
    end

    Λ, [AlgebraElem(Λ,Monomial(λGen(i))) for i=1:NLAMBDA]
end
