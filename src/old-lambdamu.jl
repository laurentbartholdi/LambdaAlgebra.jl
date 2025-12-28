
module LambdaMu

using TimerOutputs, PrettyTables, Base.Threads

export lambda_algebra

import ..LambdaAlgebra: GF, to, subscript_string, superscript_string
import ..LambdaAlgebra: degree, is_homogeneous, dimension, differential, remove_phantoms!, mark_differentials!, homology!, print_grid, print_curtis_table

################################################################
# the trigrading, in μ,λ,top
# top is the topological degree
# v+λ is the homological degree
const Degree = NamedTuple{(:μ,:λ,:top)}{NTuple{3,Int}}
Base.:+(a::Degree,b::Degree) = Degree(values(a).+values(b))
Base.:-(a::Degree,b::Degree) = Degree(values(a).-values(b))

################################################################
# Gen, λ or μ
const NGEN = 100
const NLAMBDA = NGEN÷2
const NMU = (NGEN-1)÷2

struct Gen
    x::Int16
    Gen(x) = (@assert 1≤x≤NGEN; new(x))
end

Gen(s,i) = s==:λ ? Gen(2i) : s==:μ ? Gen(2i+1) : error("Invalid generator symbol $s")
index(g::Gen) = g.x÷2
is_lambda(g::Gen) = iseven(g.x)
Base.show(io::IO, g::Gen) = print(io, is_lambda(g) ? "λ" : "μ", subscript_string(index(g)))

Base.getindex(a::AbstractArray,g::Gen...) = getindex(a,(s.x for s=g)...)
Base.setindex!(a::AbstractArray,v,g::Gen...) = setindex!(a,v,(s.x for s=g)...)
Base.convert(::Type{Gen},i::Int) = Gen(i)
Base.isless(g::Gen,h::Gen) = isless(g.x,h.x)

################################################################
# Monomial
# a vector of Gen, representing a monomial in the Lambda algebra.
# it has to be a legal expression, i.e. first v's then λ's, with
# the v's decreasing and the λ's not increasing too much.

struct Monomial <: AbstractVector{Gen}
    v::Vector{Gen}
end
Monomial() = Monomial([])
Monomial(i::Gen) = Monomial([i])
Base.length(m::Monomial) = length(m.v)
Base.size(m::Monomial) = size(m.v)
Base.:*(g::Gen,m::Monomial) = Monomial([g;m])
Base.:*(m::Monomial,n::Monomial) = Monomial(vcat(m.v,n.v))
Base.getindex(m::Monomial,i::Int) = m.v[i]
Base.getindex(m::Monomial,r::AbstractRange) = Monomial(m.v[r])
Base.setindex!(m::Monomial,g::Gen,i) = setindex!(m.v,g,i)
Base.iterate(m::Monomial,pos=1) = iterate(m.v,pos)
Base.:(==)(m::Monomial, n::Monomial) = m.v==n.v
Base.isless(m::Monomial, n::Monomial) = isless(m.v,n.v)
Base.hash(m::Monomial, h::UInt64) = hash(m.v, h)
Base.copy(m::Monomial) = Monomial(copy(m.v))

function Base.show(io::IO, m::Monomial)
    if length(m)==0
        print(io,"𝟙")
        return
    end
    i = 1
    while i≤length(m.v)
        print(io, m[i])
        exp = 1
        while i+exp≤length(m) && m[i+exp] == m[i] exp += 1 end
        exp>1 && print(io, superscript_string(exp))
        i += exp
    end
end
Base.show(io::IO, ::MIME"text/plain", m::Monomial) = show(io, m)

################################################################
# Basis
const LASTPAGE = typemax(Int)
const BASISENTRY = @NamedTuple{v::Monomial,tag::Pair{Monomial,GF{p}},phantom::Bool,page::Int} where p
const NOTAG = Monomial()
# v is the basis vector.
# page is the page at which v was tagged; LASTPAGE means still alive, 0 means phantom.
# tag is either (w,0) if v->w is a tag, or (w,c) if w->v is a tag with coefficient c.
is_alive(x::BASISENTRY) = x.page==LASTPAGE
is_phantom(x::BASISENTRY) = x.phantom
is_tagger(x::BASISENTRY) = !is_alive(x) && iszero(x.tag.second)
is_taggee(x::BASISENTRY) = !isalive(x) && !iszero(x.tag.second)

struct Basis{p} <: AbstractVector{BASISENTRY}
    data::Vector{@NamedTuple{v::Monomial,tag::Pair{Monomial,GF{p}},phantom::Bool,page::Int}}
    lookup::Dict{Monomial,Int}
    degree::Degree
end
Basis{p}(degree::Degree) where p = Basis{p}([],Dict(),degree)

function Base.show(io::IO, b::Basis)
    print(io,[x.v for x=b])
end
function Base.show(io::IO, ::MIME"text/plain", b::Basis)
    println(io,length(b),"-element basis in degree ",b.degree)
    for x=b
        print(io,"  ",x.v)
        if is_phantom(x)
            print(io, " 👻")
        elseif !is_alive(x)
            if is_tagger(x)
                print(io," → ",x.tag.first)
            else
                print(io," ← ",x.tag.second,"⋅",x.tag.first)
            end
            print(io,"[",x.page,"]")
        end
    end
end

Base.copy(b::Basis{p}) where p = Basis{p}(copy(b.data),copy(b.lookup),b.degree)
Base.iterate(b::Basis,pos=1) = iterate(b.data,pos)
Base.size(b::Basis) = size(b.data)
Base.getindex(b::Basis,i::Integer) = b.data[i]

function Base.getindex(b::Basis,m::Monomial)
    i = get(b.lookup,m,nothing)
    if i==nothing
        0,nothing
    else
        i,b.data[i]
    end
end

function Base.push!(b::Basis{p},m::Monomial,phantom::Bool) where p
    push!(b.data,(v=m,tag=(NOTAG=>zero(GF{p})),phantom=phantom,page=LASTPAGE))
    push!(b.lookup,m=>length(b.data))
    b
end

function Base.deleteat!(b::Basis, where)
    deleteat!(b.data,where)
    empty!(b.lookup)
    for i=1:length(b.data)
        push!(b.lookup,b.data[i].v=>i)
    end
    b
end

function gettag(b::Basis,m::Monomial)
    i,x = b[m]
    if i==0 || is_alive(x)
        nothing
    else
        x.tag
    end
end

function settag!(b::Basis{p},i::Int,tag::Pair{Monomial,GF{p}},page::Int=-1) where p
    b.data[i] = (v=b.data[i].v,tag=tag,phantom=b.data[i].phantom,page=(page == -1 ? b.data[i].page : page))
end

################################################################
# Algebra
struct Algebra{p}
    rules::Matrix{Union{Nothing,Vector{Tuple{GF{p},Gen,Gen}}}}
    diff::Vector{Vector{Tuple{GF{p},Gen,Gen}}}
    degree::Vector{Degree}
    basis::Dict{Degree,Basis}
    total::Ref{Int} # total degree to which we computed
    page::Ref{Int} # page number along spectral sequence
    dimension::Int
end

function Base.show(io::IO, Λ::Algebra{p}) where p
    print(io, "Λ algebra over 𝔽",subscript_string(p)," for S",superscript_string(Λ.dimension),", total degree≤",Λ.total[],", page E",superscript_string(Λ.page[]))
end
Base.show(io::IO, ::MIME"text/plain", Λ::Algebra) = show(io, Λ)

function Base.copy(Λ::Algebra{p}) where p
    Algebra{p}(Λ.rules,Λ.diff,Λ.degree,Dict(k=>copy(b) for (k,b)=Λ.basis),Ref(Λ.total[]),Ref(Λ.page[]),Λ.dimension)
end

################################################################
# AlgebraElem
struct AlgebraElem{p}
    parent::Algebra{p}
    w::Dict{Monomial,GF{p}}
end
AlgebraElem(parent::Algebra{p},w) where p = AlgebraElem{p}(parent,w)

function Base.show(io::IO, x::AlgebraElem)
    first = true
    for (k,v)=x.w
        iszero(v) && continue
        first || print(io,"+")
        isone(v) || print(io,v,"⋅")
        print(io,k)
        first = false
    end
    first && print(io,"𝟘")
end
Base.show(io::IO, ::MIME"text/plain", x::AlgebraElem) = show(io, x)

Base.iszero(x::AlgebraElem) = all(iszero,values(x.w))
Base.zero(x::AlgebraElem) = AlgebraElem(x.parent,empty(x.w))
Base.one(x::AlgebraElem{p}) where p = AlgebraElem(x.parent,Dict([]=>one(GF{p})))
Base.zero(Λ::Algebra) = AlgebraElem(Λ,Dict())
Base.one(Λ::Algebra{p}) where p = AlgebraElem(Λ,Dict([]=>one(GF{p})))

function degree(Λ::Algebra,m::Monomial)
    d = (μ=0,λ=0,top=0)
    for i=m
        d = d + Λ.degree[i]
    end
    d
end

function is_homogeneous(x::AlgebraElem)
    d = nothing    
    for k=keys(x.w)
        newd = degree(x.parent,k)
        if d==nothing
            d = newd
        elseif d≠newd
            return false
        end
    end
    return true
end

function degree(x::AlgebraElem)
    isempty(x.w) && error("0 element has no degree")
    is_homogeneous(x) || error("element is not homogeneous")
    degree(first(keys(x.w)))
end

# the first sphere dimension at which this appears
function dimension(Λ::Algebra,m::Monomial)
    isempty(m) ? 0 : m[1].x
end

dim_admissible(Λ::Algebra,m::Monomial) = dimension(Λ,m)≤Λ.dimension

function Base.:+(x::AlgebraElem{p},y::AlgebraElem{p}) where p
    @assert x.parent == y.parent
    AlgebraElem(x.parent,mergewith(+,x.w,y.w))
end
function Base.:-(x::AlgebraElem{p},y::AlgebraElem{p}) where p
    x + (-y)
end
function Base.:-(x::AlgebraElem)
    AlgebraElem(x.parent,Dict(k=>-v for (k,v)=x.w))
end
Base.:^(x::AlgebraElem,y::Integer) = Base.power_by_squaring(x,y)
Base.:*(g::GF{p},x::AlgebraElem{p}) where p = AlgebraElem(x.parent,Dict(k=>g*v for (k,v)=x.w))
Base.:*(g::Integer,x::AlgebraElem) = AlgebraElem(x.parent,Dict(k=>g*v for (k,v)=x.w))

function check(x::AlgebraElem)
    for u=keys(x.w)
        for k=1:length(u)-1
            if x.parent.rules[u[k],u[k+1]]≠nothing
                error("word $u should have been reduced")
            end
        end
    end
    return true
end

is_admissible_pair(Λ::Algebra{p},g::Gen,h::Gen) where p = p*index(g)≥index(h)+Int(is_lambda(g))

is_admissible_product(Λ::Algebra,g::Gen,m::Monomial) = isempty(m) ? is_lambda(g) : is_admissible_pair(Λ,g,m[1])

"""add_monomial!(x::AlgebraElem,u,c,from=1)

add `c*copy(u)` to `x`.
 `u` is manipulated, to avoid allocation, but is preserved on exit.
`u` is guaranteed, on entry, to be reduced at least at positions `1:2`, `2:3`, ...,`from-1:from`.
"""
function add_monomial!(x::AlgebraElem{p},u::Monomial,c::GF{p},from=1) where p
    rules = x.parent.rules
    for k=from:length(u)-1
        is_admissible_pair(x.parent,u[k],u[k+1]) && continue
        fixes = rules[u[k],u[k+1]]
        if fixes≠nothing
            uk,uknext = u[k], u[k+1]
            for s=1:length(fixes)
                u[k] = fixes[s][2]
                u[k+1] = fixes[s][3]
                newc = fixes[s][1]*c
                add_monomial!(x,u,newc,max(1,from-1))
            end
            u[k],u[k+1] = uk,uknext
            return
        end
    end
    
    dict = x.w
    stdu = getkey(dict,u,nothing) # u may be mutated later, so we need to get the internal u...
    if stdu==nothing
        dict[copy(u)] = c # ... or make a copy of it
    else
        dict[stdu] += c
    end    
    x
end

"""leading_monomial(x)

Returns `v=>c` where `c*v` is the maximal-degree term of `x`
"""
function leading_monomial(x::AlgebraElem{p}) where p
    lm = nothing
    for kv=x.w
        iszero(kv.second) && continue
        if lm == nothing || kv.first>lm.first
            lm = kv
        end
    end
    lm
end

"""tagger(Λ, m)

Searches through the bases of `Λ` to find the element n tagging `m`,
by selecting the minimal prefix of `m` such that the corresponding
suffix is in a basis; then return the prefix followed by the tagger
of the suffix.
"""
function deep_tagger(Λ::Algebra{p}, m::Monomial, deg::Degree) where p
    sign = one(GF{p})
    for i=1:length(m)-1
        deg = deg - Λ.degree[m[i]]
        if is_lambda(m[i])
            sign = -sign
        end
        γ = gettag(Λ.basis[deg],m[i+1:end])
        if γ≠nothing
            return m[1:i]*γ.first => sign*γ.second
        end
    end
end

function Base.:*(x::AlgebraElem,y::AlgebraElem)
    Λ = x.parent
    @assert Λ == y.parent
    result = zero(x)
    for (xk,xv)=x.w, (yk,yv)=y.w
        @timeit_debug to "add_monomial(*)" add_monomial!(result,xk*yk,xv*yv)
    end
    result
end

function add_differential!(result::AlgebraElem{p}, u::Monomial, sign::GF{p}) where p
    diff = result.parent.diff
    u = Gen(1)*u # prepare some space to store differential
    for s=2:length(u) # u already has an extra letter to store the differential of a letter
        ui = u[s]
        for (c,x,y)=diff[ui]
            u[s-1],u[s] = x,y
            @timeit_debug to "add_monomial(∂)" add_monomial!(result,u,sign*c)
        end
        u[s-1] = ui # put it one slot before, continue
        if is_lambda(ui)
            sign = -sign
        end
    end
    result
end
add_differential!(result::AlgebraElem{p}, m::Monomial) where p = add_differential!(result, m, one(GF{p}))

differential(Λ::Algebra{p}, m::Monomial) where p = add_differential!(zero(Λ),m)
function differential(x::AlgebraElem)
    result = zero(x)
    for lm=x.w
        add_differential!(result,lm...)
    end
    result
end

function tag_basis!(Λ::Algebra{p}, source::Basis{p}, range::Basis{p}) where p
    for (i,x)=enumerate(source) # iterate on source basis vectors
        is_alive(x) || continue
        
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

"""prebasis(Λ, degree)

Returns a basis of the `degree`-component of `Λ`, assuming all the lower-total-degree bases
have been computed. We have not yet marked the tags in it.
"""
function prebasis(Λ::Algebra{p},degree::Degree) where p
    @assert all(≥(0),degree)

    b = Basis{p}(degree)
        
    if degree.μ==degree.λ==0
        if degree.top==0
            push!(b, Monomial(),false)
        end
        return b
    end

    for g::Gen=1:NGEN
        newdegree = degree - Λ.degree[g]
        any(<(0),newdegree) && continue
        basis = Λ.basis[newdegree]
        for x=basis
            β = x.v
            is_admissible_product(Λ,g,β) || continue
            if !is_alive(x)
                γ = x.tag.first
                is_admissible_product(Λ,g,γ) && continue
            end
            
            γ = g*β
            push!(b,γ,!dim_admissible(Λ,γ))
        end
    end
    @assert issorted([x.v for x=b])
    b
end

"""cache_basis(Λ, total)

Caches the basis of `Λ` up to total degree `total`, and marks its tags.
"""
function cache_basis!(Λ::Algebra,total::Int)
    while Λ.total[]<total
        Λ.total[] += 1

        for hom=0:Λ.total[]
            top = Λ.total[]-hom
            for μ=0:hom
                degree = (μ=μ,λ=hom-μ,top=top)
                Λ.basis[degree] = prebasis(Λ,degree)
            end
        end

        # we can do this in parallel, since all μ's operate independently
        #@threads # speedup seems to be only 50%
        for μ=0:Λ.total[]-1
            for λ=0:Λ.total[]-μ-1
                top = Λ.total[]-λ-μ
                source = (μ=μ,λ=λ,top=top)
                range = (μ=μ,λ=λ+1,top=top-1)
                tag_basis!(Λ,Λ.basis[source],Λ.basis[range])
            end
        end
    end
end

function remove_phantoms!(Λ::Algebra{p},b::Basis{p}) where p
    b.degree==(μ=2,λ=2,top=14) && @info "rem" b.degree [is_phantom(x) || !dim_admissible(Λ,x.v) for x=b]
    deleteat!(b,[is_phantom(x) || !dim_admissible(Λ,x.v) for x=b])
    b.degree==(μ=2,λ=2,top=14) && @info "rem2" b.data    
    # now manually erase tags that point to nowhere. could be optimized.
    for i=1:length(b)
        if !is_alive(b[i]) && !dim_admissible(Λ,b[i].tag.first)
            b.data[i] = (v=b[i].v,tag=(NOTAG=>zero(GF{p})),phantom=false,page=LASTPAGE)
        end
    end
end

"""remove_phantoms!(Λ,[b::Basis])

Remove all basis entries that have dimension too large
"""
function remove_phantoms!(Λ::Algebra)
    foreach(values(Λ.basis)) do b
        remove_phantoms!(Λ,b)
    end
    Λ
end

"""homology(Λ)

Removes all source and range monomials from the differential.

Returns the number of removed differentials.
"""
function homology!(Λ::Algebra)
    count = 0

    foreach(values(Λ.basis)) do b
        cancelled = [!is_alive(x) for x=b]
        deleteat!(b,cancelled)
        count += Base.count(cancelled)
    end
    count
end

function mark_differentials!(source::Vector{Basis{p}},range::Vector{Basis{p}}) where p
    M = length(source)
    @assert M==length(range)

    # for all vectors in all source bases, in increasing order:
    # either give them a tag to range basis vector, & back, or mark them as untagged.
    # we steal the bit "phantom" to mark a source basis vector that we already consulted

    @assert all(b->!any(b.phantom),source) # we'll need these to keep track of untagged vectors

    count = 0
    
    for i=1:M
        for k=1:length(source[i])
            v = source[i][k]
            source[i].phantom[k] && continue # already consulted, should not be tagged
            source[i].tagger[k]==source[i].taggee[k]==nothing || continue # already tagged
            taggees = [range[i][l] for l=1:length(range[i]) if range[i].tagger[l]==range[i].taggee[l]==nothing]
            for j=i+1:M
                v∈source[j] || setdiff!(taggees,range[j])
            end
            if isempty(taggees)
                for j=i:M
                    l = source[j][v]
                    if l≠nothing
                        source[j].phantom[l] = true
                    end
                end
            else
                w = first(taggees)
                for j=i:M
                    l = source[j][v]
                    m = range[j][w]
                    if l≠nothing && m≠nothing
                        count += 1
                        source[j].taggee[l] = w
                        range[j].tagger[m] = v=>one(GF{p})
                    end
                end
            end
        end
    end

    # reset the "phantom" bits
    for i=1:M
        fill!(source[i].phantom,false)
    end

    count
end

"""mark_differentials!(Λs::Vector{Algebra},Λ₊,dμ,dλ)

Add tags for higher differentials in each `Λ[n]`, with tridegree `(dμ,dλ,-1)`, that are
compatible with the suspension maps `Λ[n] → Λ[n+1]`.

Returns the number of new differentials.
"""
function mark_differentials!(Λ::Vector{Algebra{p}}, d::Degree) where p
    total = Λ[1].total[]
    @assert all(L->L.total[]==total,Λ)
    @assert d.top == -1
    @assert d.λ ≥ 1
    @assert d.μ ≥ 0
    
    count = 0

    for hom=total+1-d.μ-d.λ:-1:1
        for top=hom:total+1-hom-d.μ-d.λ
            for μ=0:hom
                sourcedegree = (μ=μ,λ=hom-μ,top=top)
                rangedegree = sourcedegree + d
                count += mark_differentials!([L.basis[sourcedegree] for L=Λ],[L.basis[rangedegree] for L=Λ])
            end
        end
    end

    count
end
mark_differentials!(Λ::Vector{Algebra{p}}, μ::Int, λ::Int) where p = mark_differentials!(Λ,(μ=μ,λ=λ,top=-1))

function propagate_tags!(Λ::Algebra, source::Basis, range::Basis)
    sourcedel = Int[]
    rangedel = Int[]
    for i=1:length(source) # iterate on source basis vectors
        source.taggee[i]==nothing && continue
        β = source[i]
        γ = source.taggee[i]
        j = range[γ]
        sourcedegree = source.degree
        rangedegree = range.degree
        for s=1:min(length(β),length(γ))
            β[s]==γ[s] || break
            sourcedegree -= Λ.degree[β[s]]
            rangedegree -= Λ.degree[β[s]]
            k = Λ.basis[sourcedegree][β[s+1:end]]
            k==nothing && continue
            Λ.basis[sourcedegree].taggee[k]==nothing && continue
            if Λ.basis[sourcedegree].taggee[k]==γ[s+1:end]
                push!(sourcedel,i)
                push!(rangedel,j)
                break
            end
        end
    end
    deleteat!(source,sourcedel)
    deleteat!(range,sort(rangedel))
end

"""propagate_tags!(Λ,[source::Basis,range::Basis])

removes all pairs (tagger,taggee) that are consequences of shorter pairs
"""
function propagate_tags!(Λ::Algebra)
    for total=0:Λ.total[]
        for hom=0:total-1
            top = total-hom
            for μ=0:hom
                source = (μ=μ,λ=hom-μ,top=top)
                range = (μ=μ,λ=hom-μ+1,top=top-1)
                propagate_tags!(Λ,Λ.basis[source],Λ.basis[range])
            end
        end
    end
    Λ
end

function Base.getindex(Λ::Algebra,degree::Degree)
    cache_basis!(Λ,degree.μ+degree.λ+degree.top)
    Λ.basis[degree]
end

function Base.getindex(Λ::Algebra,i::Int,j::Int,top::Int)
    cache_basis!(Λ,i+j+top)
    Λ[(μ=i,λ=j,top=top)]
end

function Base.getindex(Λ::Algebra,::Colon,j::Int,top::Int)
    vcat((Λ.basis[(μ=μ,λ=j,top=top)].vectors for μ=0:Λ.total[]-j-top)...)
end

function Base.getindex(Λ::Algebra,ij::Int,top::Int)
    cache_basis!(Λ,ij+top)
    vcat((Λ.basis[(μ=μ,λ=ij-μ,top=top)].vectors for μ=0:ij)...)
end

function Base.getindex(Λ::Algebra,top::Int)
    vcat((Λ.basis[(μ=i,λ=j,top=top)].vectors for i=0:top for j=0:top-i)...)
end

"""lambda_algebra(p=3)

Construct the lambda algebra in odd characteristic `p`
"""
function lambda_algebra(;p=3, dimension=99, top_degree=-1)
    @assert isodd(p)
    K = GF{p}

    A(k,j) = K((-1)^(j+1) * binomial(big(p-1)*(k-j)-1,j))
    B(k,j) = K((-1)^j * binomial(big(p-1)*(k-j),j))
    
    rules = Matrix{Union{Nothing,Vector{Tuple{K,Gen,Gen}}}}(undef,NGEN,NGEN)
    fill!(rules,nothing)
    for i=1:NLAMBDA, k=0:NLAMBDA-p*i # λᵢλⱼ
        rule = []
        for j=0:k-1
            iszero(A(k,j)) || push!(rule,(A(k,j),Gen(:λ,i+k-j),Gen(:λ,p*i+j)))
        end
        rules[Gen(:λ,i),Gen(:λ,p*i+k)] = rule
    end
    for i=1:NLAMBDA, k=0:NMU-p*i # λᵢμⱼ
        rule = []
        for j=0:k-1
            iszero(A(k,j)) || push!(rule,(A(k,j),Gen(:λ,i+k-j),Gen(:μ,p*i+j)))
        end
        for j=0:k
            iszero(B(k,j)) || push!(rule,(B(k,j),Gen(:μ,i+k-j),Gen(:λ,p*i+j)))
        end
        rules[Gen(:λ,i),Gen(:μ,p*i+k)] = rule
    end
    for i=0:NMU, k=0:NLAMBDA-p*i-1 # μᵢλⱼ
        rule = []
        for j=0:k-1
            iszero(A(k,j)) || push!(rule,(A(k,j),Gen(:μ,i+k-j),Gen(:λ,p*i+j+1)))
        end
        rules[Gen(:μ,i),Gen(:λ,p*i+k+1)] = rule
    end
    for i=0:NMU, k=0:NMU-p*i-1 # μᵢμⱼ
        rule = []
        for j=0:k-1
            iszero(A(k,j)) || push!(rule,(A(k,j),Gen(:μ,i+k-j),Gen(:μ,p*i+j+1)))
        end
        rules[Gen(:μ,i),Gen(:μ,p*i+k+1)] = rule
    end
    
    diff = Vector{Vector{Tuple{K,Gen,Gen}}}(undef,NGEN)
    for k=1:NLAMBDA
        diff[Gen(:λ,k)] = [(A(k,j),Gen(:λ,k-j),Gen(:λ,j)) for j=1:k-1 if !iszero(A(k,j))]
    end
    for k=0:NMU
        diff[Gen(:μ,k)] = [[(A(k,j),Gen(:λ,k-j),Gen(:μ,j)) for j=0:k-1 if !iszero(A(k,j))];
                           [(B(k,j),Gen(:μ,k-j),Gen(:λ,j)) for j=1:k if !iszero(B(k,j))]]
    end
    
    degrees = Vector{Degree}(undef,NGEN)
    for i=1:NLAMBDA
        degrees[Gen(:λ,i)] = (μ=0,λ=1,top=(2p-2)*i-1)
    end
    for j=0:NMU
        degrees[Gen(:μ,j)] = (μ=1,λ=0,top=(2p-2)*j)
    end
                                            
    Λ = Algebra{p}(rules,diff,degrees,Dict(),Ref(-1),Ref(1),dimension)
    cache_basis!(Λ,top_degree)
    
    # sanity check
    for g::Gen=1:NGEN, h::Gen=1:NGEN
        @assert is_admissible_pair(Λ,g,h) == (rules[g,h]==nothing)
    end

    Λ, [AlgebraElem(Λ,Dict(Monomial(Gen(:λ,i))=>one(K))) for i=1:NLAMBDA],[AlgebraElem(Λ,Dict(Monomial(Gen(:μ,i))=>one(K))) for i=0:NMU]
end

function print_curtis_table(io::IO, Λ::Algebra, totaldegree::Int)
    curtis_string(g) = string(is_lambda(g) ? index(g) : -index(g))
    cache_basis!(Λ,totaldegree)
    for top=1:totaldegree
        for hom=1:totaldegree-top
            first = true
            for μ=hom:-1:0
                basis = Λ.basis[(μ=μ,λ=hom-μ,top=top)]
                for i=length(basis):-1:1
                    β = basis[i]
                    if basis.taggee[i]==nothing
                        if first
                            print(io, top,",",hom)
                            first = false
                        end
                        print(io," 1(",join((curtis_string(g) for g=β.v)," "),")")
                        if basis.tagger[i]≠nothing
                            γ = basis.tagger[i]
                            print(io, "/",γ.second,"(", join((curtis_string(g) for g=γ.first)," "),")")
                        end
                    end
                end
            end
            first || println(io)
        end
    end
end
print_curtis_table(Λ::Algebra, totaldegree::Int) = print_curtis_table(stdout,Λ,totaldegree)

function print_grid(io::IO, Λ::Algebra{p}; top_degree = nothing, names = nothing, width = nothing, kwargs...) where p
    if isa(top_degree,Int)
        cache_basis!(Λ,top_degree + top_degree÷(2p-2) + 1)
    else
        top_degree = Λ.total[]
    end
    if names==nothing
        names = top_degree ≤ 30
    end
    if width==nothing
        width = names ? 150 : 60
    end
    
    backend = get(kwargs,:backend,:text)     

    data = Matrix{String}(undef,top_degree,top_degree)
    maxhom = 1
    fill!(data,"")
    arrows = Pair{NTuple{3,Int},NTuple{3,Int}}[]
    for top=3:top_degree
        for hom=1:top
            if hom+top > Λ.total[]
                data[top,hom] = "?"
                continue
            end
            entries = String[]
            for μ=0:hom
                basis = Λ.basis[(μ=μ,λ=hom-μ,top=top)]
                i = 0
                for x=basis
                    is_phantom(x) && continue
                    i += 1
                    if is_tagger(x) && dim_admissible(Λ,x.tag.first) && backend==:tikz
                        rangedegree = degree(Λ,x.tag.first)
                        j = 0
                        for y=Λ.basis[rangedegree] # find position in range
                            is_phantom(y) && continue
                            j += 1
                            y.v==x.tag.first && break
                        end
                        push!(arrows,(top,hom,i)=>(top-1,rangedegree.λ+rangedegree.μ,j))
                    end
                    if is_tagger(x) && dim_admissible(Λ,x.tag.first) && backend==:html
                        rangedegree = degree(Λ,x.tag.first)
                        shift = rangedegree.λ+rangedegree.μ-(basis.degree.λ+basis.degree.μ)
                        x1,x2 = 5, 5+width*shift
                        y1,y2 = 12,28
                        arrow = """<svg width="$x2" style="position: absolute;
    width: 0;
    height: 0;
    top: 0;
    left: 0;">
  <defs>
    <marker id="redhead" orient="auto" markerWidth="6" markerHeight="4" refX="6" refY="2" orient="auto" markerUnits="strokeWidth">
      <path d="M0,0 L0,4 L6,2 Z" fill="red"/>
    </marker>
  </defs>
  <line x1="$x1" y1="$y1" x2="$x2" y2="$y2" style="stroke:rgb(255,0,0);stroke-width:1.5" marker-end="url(#redhead)"/>
</svg>"""
                    else
                        arrow = ""
                    end
                    if isempty(arrow)                        
                        if names
                            s = string(x.v)
                        else
                            s = string(μ)
                        end
                        if backend==:html
                            push!(entries, """<div style="background-color:#d5f4e6;">$s</div>""")
                        else
                            push!(entries, s)
                        end
                    else
                        if names
                            push!(entries, arrow*string(x.v))
                        else
                            push!(entries, arrow*string(μ))
                        end
                    end
                end
            end
            if !isempty(entries)
                maxhom = max(hom,maxhom)
                data[top,hom] = join(entries,",")
            end
        end
    end
    @info entries
    if backend==:latex
        kwargs = (highlighters = [LatexHighlighter((_, i, _)->iszero((top_degree+1-i)%5), ["textbf"])],
                  kwargs...)
    elseif backend==:html
        kwargs = (maximum_column_width = string(width,"px"),
                  allow_html_in_cells = true,
                  style = HtmlTableStyle(first_line_column_label = ["width" => string(width,"px")]),
                  kwargs...)
    elseif backend==:markdown
        kwargs = (highlighters = [MarkdownHighlighter((_, i, _)->iszero((top_degree+1-i)%5), MarkdownStyle(bold=true))],
                  kwargs...)
    else 
        kwargs = (display_size = (-1,-1),
                  
                  highlighters = [TextHighlighter((_, i, _)->iszero((top_degree+1-i)%5), crayon"fg:black bold bg:light_gray")],
                  #table_format = TextTableFormat(horizontal_lines_at_data_rows=collect(mod1(top_degree,5):5:top_degree)),
                  kwargs...)
    end

    if backend==:tikz
        println(io,"""\\documentclass{standalone}
\\usepackage{tikz}
\\usetikzlibrary{matrix,positioning,calc}
\\usepackage{unicode-math}
\\begin{document}
\\tikzset{toprule/.style={%
        execute at end cell={%
            \\draw [line cap=rect,#1] (\\tikzmatrixname-\\the\\pgfmatrixcurrentrow-\\the\\pgfmatrixcurrentcolumn.north west) -- (\\tikzmatrixname-\\the\\pgfmatrixcurrentrow-\\the\\pgfmatrixcurrentcolumn.north east);%
        }
    },
    bottomrule/.style={%
        execute at end cell={%
            \\draw [line cap=rect,#1] (\\tikzmatrixname-\\the\\pgfmatrixcurrentrow-\\the\\pgfmatrixcurrentcolumn.south west) -- (\\tikzmatrixname-\\the\\pgfmatrixcurrentrow-\\the\\pgfmatrixcurrentcolumn.south east);%
        }
    }
}
\\begin{tikzpicture}[node distance=0.5em,>=stealth]
\\matrix[row sep=2ex,column sep=1em] {""")
        for i=1:maxhom
            print(io," & \\node{\\textbf{",i,"}};")
        end
        println(io, "\\\\")
        for i=top_degree:-1:3
            print(io,"\\node[minimum height=3ex,rectangle] (c",i,"-begin) {\\textbf{",i,"}};")
            for j=1:maxhom
                s = split(data[i,j],",")
                print(io,"&")
                placement = ""
                count = 1
                for a=s
                    print(io," \\node[inner sep=0pt,outer sep=0pt,",placement,"] (c",i,"-",j,"-",count,") {\\(",a,"\\)};")
                    placement = "right=of c$i-$j-$count"
                    count += 1
                end
                println(io,"\\node[",placement,"] (c",i,"-",j,"-end) {};")
            end
            println(io,"\\\\")
        end
        println(io,"};")
        for (a,b)=arrows
            println(io,"\\draw[red,->] (c",a[1],"-",a[2],"-",a[3],".south) -- (c",b[1],"-",b[2],"-",b[3],".north);")
        end
        for i=mod1(top_degree,5):5:top_degree
            println(io,"\\draw (c$i-begin.north west) -- (c$i-$maxhom-end.north east |- c$i-begin.north west);")
        end
        println(io,"""\\end{tikzpicture}
\\end{document}""")
        return
    end
    
    pretty_table(io, data[top_degree:-1:3,1:maxhom];
                 column_labels=1:maxhom, row_labels=top_degree:-1:3,
                 title = "E² page of $Λ",
                 kwargs...)
end
print_grid(Λ::Algebra; kwargs...) = print_grid(stdout, Λ; kwargs...)

function pdf_grid(name::AbstractString, Λ::Algebra{p}; top_degree = nothing, names = nothing) where p
    mktempdir() do dir
        cd(dir) do
            tex = tempname(dir,suffix=".tex")
            open(tex,"w") do f
                print_grid(f, Λ, top_degree=top_degree, names=names, backend=:tikz)
            end
            out = IOBuffer()
            if !success(pipeline(`xelatex $tex`,stdout=out))
                @error String(take!(out))
            end
            mv(tex[1:end-3]*"pdf",name,force=true)
        end
    end
end

end

                                     
