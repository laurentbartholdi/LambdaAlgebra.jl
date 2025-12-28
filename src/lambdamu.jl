const LAMBDAMU = (:λ,:μ)

index(g::Gen{p,LAMBDAMU}) where p = g.x÷2
kind(g::Gen{p,LAMBDAMU}) where p = iseven(g.x) ? 'λ' : 'μ'

is_admissible_pair(g::Gen{p,LAMBDAMU},h::Gen{p,LAMBDAMU}) where {p} = p*index(g)≥index(h)+Int(is_lambda(g))
is_admissible_product(g::Gen{p,LAMBDAMU},m::Monomial{p,LAMBDAMU}) where {p} = isempty(m) ? is_lambda(g) : is_admissible_pair(g,m[1])

# the first sphere dimension at which this appears
dimension(Λ::Algebra{p,LAMBDAMU},m::Monomial{p,LAMBDAMU}) where {p} = isempty(m) ? 0 : m[1].x

"""lambda_algebra(p=3)

Construct the lambda algebra in odd characteristic `p`
"""
function lambda_algebra(;p=3, top_degree=-1, dimension=STABLE_DIMENSION)
    p==2 && return lambda2_algebra(;top_degree,dimension)
    
    @assert isodd(p) # actually test isprime
    K = GF{p}

    A(k,j) = K((-1)^(j+1) * binomial(big(p-1)*(k-j)-1,j))
    B(k,j) = K((-1)^j * binomial(big(p-1)*(k-j),j))

    GEN = Gen{p,LAMBDAMU}
    NLAMBDA = NGEN÷2
    NMU = (NGEN-1)÷2
    
    λGen(i) = GEN(2i)
    μGen(i) = GEN(2i+1)
    
    rules = Matrix{Union{Nothing,Vector{Tuple{K,GEN,GEN}}}}(undef,NGEN,NGEN)
    fill!(rules,nothing)
    for i=1:NLAMBDA, k=0:NLAMBDA-p*i # λᵢλⱼ
        rule = []
        for j=0:k-1
            iszero(A(k,j)) || push!(rule,(A(k,j),λGen(i+k-j),λGen(p*i+j)))
        end
        rules[λGen(i),λGen(p*i+k)] = rule
    end
    for i=1:NLAMBDA, k=0:NMU-p*i # λᵢμⱼ
        rule = []
        for j=0:k-1
            iszero(A(k,j)) || push!(rule,(A(k,j),λGen(i+k-j),μGen(p*i+j)))
        end
        for j=0:k
            iszero(B(k,j)) || push!(rule,(B(k,j),μGen(i+k-j),λGen(p*i+j)))
        end
        rules[λGen(i),μGen(p*i+k)] = rule
    end
    for i=0:NMU, k=0:NLAMBDA-p*i-1 # μᵢλⱼ
        rule = []
        for j=0:k-1
            iszero(A(k,j)) || push!(rule,(A(k,j),μGen(i+k-j),λGen(p*i+j+1)))
        end
        rules[μGen(i),λGen(p*i+k+1)] = rule
    end
    for i=0:NMU, k=0:NMU-p*i-1 # μᵢμⱼ
        rule = []
        for j=0:k-1
            iszero(A(k,j)) || push!(rule,(A(k,j),μGen(i+k-j),μGen(p*i+j+1)))
        end
        rules[μGen(i),μGen(p*i+k+1)] = rule
    end
    
    diff = Vector{Vector{Tuple{K,GEN,GEN}}}(undef,NGEN)
    for k=1:NLAMBDA
        diff[λGen(k)] = [(A(k,j),λGen(k-j),λGen(j)) for j=1:k-1 if !iszero(A(k,j))]
    end
    for k=0:NMU
        diff[μGen(k)] = [[(A(k,j),λGen(k-j),μGen(j)) for j=0:k-1 if !iszero(A(k,j))];
                           [(B(k,j),μGen(k-j),λGen(j)) for j=1:k if !iszero(B(k,j))]]
    end
    
    degrees = Vector{Degree}(undef,NGEN)
    for i=1:NLAMBDA
        degrees[λGen(i)] = (μ=0,λ=1,top=(2p-2)*i-1)
    end
    for j=0:NMU
        degrees[μGen(j)] = (μ=1,λ=0,top=(2p-2)*j)
    end
                                            
    Λ = Algebra{p,LAMBDAMU}(rules,diff,degrees,Dict(),Ref(-1),Ref(1),Ref(dimension))
    # sanity check
    for g::GEN=1:NGEN, h::GEN=1:NGEN
        @assert is_admissible_pair(g,h) == (rules[g,h]==nothing)
    end

    cache_basis!(Λ,top_degree)
    if dimension≠STABLE_DIMENSION
#        truncate!(Λ,abs(dimension))
    end

    Λ, [AlgebraElem(Λ,SortedDict(Monomial(λGen(i))=>one(K))) for i=1:NLAMBDA],[AlgebraElem(Λ,SortedDict(Monomial(μGen(i))=>one(K))) for i=0:NMU]
end
