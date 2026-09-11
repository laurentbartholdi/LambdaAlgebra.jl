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

# overload truncate, we have to add back some implied terms
function ___truncate!(Λ::Algebra{p,LAMBDAV},max_dimension::Int) where {p}
    @assert isodd(max_dimension)

    for total=0:Λ.total[]
        for hom=0:total
            top = total-hom
            for μ=0:hom
                b = Λ.basis[(μ=μ,λ=hom-μ,top=top)]
                todelete = Int[]
                for (i,x) = enumerate(b)
                    if is_taggee(x) && dimension(Λ,x.tag.first)>max_dimension # erase its tag
                        settag!(b,i,b.notag,LASTPAGE)
                    end
                    if dimension(Λ,x.v)>max_dimension
                        if is_tagger(x)
                            for g::Gen{p,LAMBDAV}=NGEN-NV:NGEN # only try to add v's
                                is_admissible_product(g,x.v) || break # too large v, stop
                                is_admissible_product(g,x.tag.first) || break # too large v, stop
                                @info "pushing $g to $x"
                                γ = g*x.v
                                δ = g*x.tag.first
                                sourcedegree = degree(Λ,γ)
                                rangedegree = degree(Λ,δ)
                                haskey(Λ.basis,sourcedegree) || continue # does not exist in basis
                                haskey(Λ.basis,rangedegree) || continue
                                @info "push!" (Λ.basis[sourcedegree],γ,δ)
                                push!(Λ.basis[sourcedegree],γ;tag=δ,page=Λ.page[])
                                @info "push!" (Λ.basis[rangedegree],δ,γ)
                                push!(Λ.basis[rangedegree],δ;tag=γ=>one(GF{p}),page=Λ.page[])
                                dimension(Λ,γ)≤max_dimension && break # this arrow will stay
                            end
                        end
                        push!(todelete,i)
                    end
                end
                deleteat!(b,todelete)
            end
        end
    end
    Λ
end

# overload for periodic algebra: only add dimension-OK tags
function tag_basis!(Λ::Algebra{p,LAMBDAV}, source::Basis{p,LAMBDAV}, range::Basis{p,LAMBDAV}) where {p}
    for (i,x)=enumerate(source) # iterate on source basis vectors
        is_alive(x) || continue
        dimension(Λ,x.v)≤Λ.dimension[] || continue
        
        @timeit_debug to "tag_monomial" begin
            @timeit_debug to "differential" δ = differential(Λ,x.v)
            while !iszero(δ)
                m,c = leading_monomial(δ)
                j,y = range[m]
                if j==0
                    @timeit_debug to "tagger" γ = deep_tagger(Λ,m,range.degree)
                else
                    if is_alive(y)
                        settag!(source,i,y.v=>zero(GF{p}),Λ.page[])
                        settag!(range,j,x.v=>inv(c),Λ.page[])
                        break
                    end
                    γ = y.tag
                end
                @timeit_debug to "add_differential" add_differential!(δ,γ.first,(-c)*γ.second)
            end
        end
    end
end

# overload for periodic algebra: add back some of the missing tags
function prebasis(Λ::Algebra{p,LAMBDAV},degree::Degree) where {p}
    @assert all(≥(0),degree)

    b = Basis{p,LAMBDAV}(degree)
        
    if degree.μ==degree.λ==0
        if degree.top==0
            push!(b, Monomial{p,LAMBDAV}())
        end
        return b
    end

    for g::Gen{p,LAMBDAV}=1:NGEN
        newdegree = degree - Λ.degree[g]
        any(<(0),newdegree) && continue
        basis = Λ.basis[newdegree]
        for x=basis
            β = x.v
            is_admissible_product(g,β) || continue
            γ = g*β
            if Λ.dimension[]==STABLE_DIMENSION && !is_alive(x)
                δ = x.tag.first
                is_admissible_product(g,δ) && is_lambda(g) && continue
                is_admissible_product(g,δ) && dimension(Λ,β)≤Λ.dimension[] && dimension(Λ,δ)≤Λ.dimension[] && continue
            end            
            push!(b,γ)
        end
    end
#    sort!(b.data,by=x->(dimension(Λ,x.v)>Λ.dimension[],x.v))
#    empty!(b.lookup)
#    for (i,x)=enumerate(b.data) push!(b.lookup,x.v=>i) end
    @assert issorted(b)
    b
end

"""periodic_algebra(p=3)

Construct the lambda algebra in odd characteristic `p`
"""
function periodic_algebra(;p=3, top_degree=-1, μ_degree=-1, dimension=STABLE_DIMENSION)
    @assert isodd(p)
    K = GF{p}

    NLAMBDA = NGEN-NV-1
    
    A(k,j) = K((-1)^(j+1) * binomial(big(p-1)*(k-j)-1,j))

    GEN = Gen{p,LAMBDAV}
    λGen(i) = GEN(i)
    vGen(i) = GEN(i+NGEN-NV)
    
    rules = Matrix{Union{Nothing,Vector{Tuple{K,GEN,GEN}}}}(undef,NGEN,NGEN)
    fill!(rules,nothing)
    fill!(rules,nothing)
    for i=0:NV, j=i+1:NV # all v's increasing
        rules[vGen(j),vGen(i)] = [(one(K),vGen(i),vGen(j))]
    end
    for i=1:NLAMBDA, j=0:NV # move v's before λ's
        if j==0 || i+p^(j-1) > NLAMBDA
            rules[λGen(i),vGen(j)] = [(one(K),vGen(j),λGen(i))]
        else
            rules[λGen(i),vGen(j)] = [(one(K),vGen(j),λGen(i)),(one(K),vGen(j-1),λGen(i+p^(j-1)))]
        end
    end
    for i=1:NLAMBDA, k=0:NLAMBDA-p*i # λ's don't increase too much
        rule = []
        for j=0:k-1
            iszero(A(k,j)) || push!(rule,(A(k,j),λGen(i+k-j),λGen(p*i+j)))
        end
        rules[λGen(i),λGen(p*i+k)] = rule
    end
    
    diff = Vector{Vector{Tuple{K,Gen,Gen}}}(undef,NGEN)
    for k=1:NLAMBDA
        diff[λGen(k)] = [(A(k,j),λGen(k-j),λGen(j)) for j=1:k-1 if !iszero(A(k,j))]
    end
    diff[vGen(0)] = [] # d(v₀) = 0
    for j=1:NV
        diff[vGen(j)] = [(one(K),vGen(j-1),λGen(p^(j-1)))]
    end
    
    degrees = Vector{Degree}(undef,NGEN)
    for i=1:NLAMBDA
        degrees[λGen(i)] = (μ=0,λ=1,top=(2p-2)*i-1)
    end
    for j=0:NV
        degrees[vGen(j)] = (μ=1,λ=0,top=2*p^j-2)
    end
                                            
    Λ = Algebra{p,LAMBDAV}(rules,diff,degrees,dimension)
    # sanity check
    for g::GEN=1:NGEN, h::GEN=1:NGEN
        @assert is_admissible_pair(g,h) == (rules[g,h]==nothing)
    end

    cache_basis!(Λ,top_degree,μ_degree)
    if dimension≠STABLE_DIMENSION
        truncate!(Λ,dimension)
    end
    
    Λ, [AlgebraElem(Λ,Monomial(λGen(i))) for i=1:NLAMBDA],[AlgebraElem(Λ,Monomial(vGen(i))) for i=0:NV]
end
