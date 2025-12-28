################################################################
# GF(p)
struct GF{p}
    x::UInt8
    GF{p}(x::Integer) where p = new(mod(x,p))
end
Base.Int(x::GF) = x.x
Base.show(io::IO, ::MIME"text/plain", x::GF) = print(io, x.x)
Base.show(io::IO, x::GF) = print(io, x.x)
Base.:+(x::GF{p},y::GF{p}) where p = GF{p}(x.x+y.x)
Base.:-(x::GF{p},y::GF{p}) where p = GF{p}(p+x.x-y.x)
Base.:-(x::GF{p}) where p = GF{p}(p-x.x)
Base.:*(x::GF{p},y::GF{p}) where p = GF{p}(UInt16(x.x)*y.x)
Base.:*(x::Integer,y::GF{p}) where p = GF{p}(x*y.x)
Base.inv(x::GF{p}) where p = GF{p}(gcdx(Int(x.x),p)[2])
Base.one(::Type{GF{p}}) where p = GF{p}(1)
Base.zero(::Type{GF{p}}) where p = GF{p}(0)
Base.one(::GF{p}) where p = GF{p}(1)
Base.zero(::GF{p}) where p = GF{p}(0)

################################################################
# the trigrading, in μ,λ,top
# top is the topological degree
# v+λ is the homological degree
const Degree = NamedTuple{(:μ,:λ,:top)}{NTuple{3,Int}}
Base.:+(a::Degree,b::Degree) = Degree(values(a).+values(b))
Base.:-(a::Degree,b::Degree) = Degree(values(a).-values(b))

################################################################
# Gen, λ or μ or v
const NGEN = 100

struct Gen{p,Names}
    x::Int16
    Gen{p,Names}(x) where {p,Names} = (@assert 1≤x≤NGEN; new{p,Names}(x))
end

function kind end
function index end
is_lambda(g::Gen) = kind(g)=='λ'
        
Base.show(io::IO, g::Gen{p,Names}) where {p,Names} = print(io, kind(g), subscript_string(index(g)))

Base.getindex(a::AbstractArray,g::Gen,gs::Gen...) = getindex(a,g.x,(s.x for s=gs)...)
Base.setindex!(a::AbstractArray,v,g::Gen,gs::Gen...) = setindex!(a,v,g.x,(s.x for s=gs)...)
Base.convert(T::Type{Gen{p,Names}},i::Int) where {p,Names} = T(i)
Base.isless(g::Gen,h::Gen) = isless(g.x,h.x)

################################################################
# Monomial
# a vector of Gen, representing a monomial in the Lambda algebra.
# it has to be a legal expression, i.e. first v's then λ's, with
# the v's decreasing and the λ's not increasing too much.

struct Monomial{p,Names} <: AbstractVector{Gen{p,Names}}
    v::Vector{Gen{p,Names}}
end
Monomial{p,Names}() where {p,Names} = Monomial{p,Names}([])
Monomial(i::Gen{p,Names}) where {p,Names} = Monomial([i])
Base.length(m::Monomial) = length(m.v)
Base.size(m::Monomial) = size(m.v)
Base.:*(g::Gen{p,Names},m::Monomial{p,Names}) where {p,Names} = Monomial([g;m])
Base.:*(m::Monomial{p,Names},n::Monomial{p,Names}) where {p,Names} = Monomial(vcat(m.v,n.v))
Base.getindex(m::Monomial,i::Int) = m.v[i]
Base.getindex(m::Monomial,r::AbstractRange) = Monomial(m.v[r])
Base.setindex!(m::Monomial{p,Names},g::Gen{p,Names},i) where {p,Names} = setindex!(m.v,g,i)
Base.iterate(m::Monomial,pos=1) = iterate(m.v,pos)
Base.:(==)(m::Monomial{p,Names}, n::Monomial{p,Names}) where {p,Names} = m.v==n.v
Base.isless(m::Monomial{p,Names}, n::Monomial{p,Names}) where {p,Names} = isless(m.v,n.v)
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
const BASISENTRY = @NamedTuple{v::Monomial{p,Names},tag::Pair{Monomial{p,Names},GF{p}},page::Int} where {p,Names}

# v is the basis vector.
# page is the page at which v was tagged; LASTPAGE means still alive
# tag is either (w,0) if v->w is a tag, or (w,c) if w->v is a tag with coefficient c.
is_alive(x::BASISENTRY) = x.page==LASTPAGE
is_tagger(x::BASISENTRY) = !is_alive(x) && iszero(x.tag.second)
is_taggee(x::BASISENTRY) = !is_alive(x) && !iszero(x.tag.second)

struct Basis{p,Names} <: AbstractVector{BASISENTRY}
    data::Vector{BASISENTRY}
    lookup::Dict{Monomial{p,Names},Int}
    notag::Pair{Monomial{p,Names},GF{p}} # store it once
    degree::Degree
end
Basis{p,Names}(degree::Degree) where {p,Names} = Basis{p,Names}([],Dict(),Monomial{p,Names}()=>zero(GF{p}),degree)

function Base.show(io::IO, b::Basis{p,Names}) where {p,Names}
    print(io,Monomial{p,Names}[x.v for x=b])
end
function Base.show(io::IO, ::MIME"text/plain", b::Basis)
    println(io,length(b),"-element basis in degree ",b.degree)
    for x=b
        print(io,"  ",x.v)
        if !is_alive(x)
            if is_tagger(x)
                print(io," → ",x.tag.first)
            else
                print(io," ← ",x.tag.second,"⋅",x.tag.first)
            end
            print(io,"[",x.page,"]")
        end
    end
end

Base.copy(b::Basis{p,Names}) where {p,Names} = Basis{p,Names}(copy(b.data),copy(b.lookup),b.notag,b.degree)
Base.iterate(b::Basis,pos=1) = iterate(b.data,pos)
Base.size(b::Basis) = size(b.data)
Base.getindex(b::Basis,i::Integer) = b.data[i]

function Base.getindex(b::Basis{p,Names},m::Monomial{p,Names}) where {p,Names}
    i = get(b.lookup,m,nothing)
    if i==nothing
        0,nothing
    else
        i,b.data[i]
    end
end

function Base.in(m::Monomial{p,Names},b::Basis{p,Names}) where {p,Names}
    haskey(b.lookup,m)
end

function Base.push!(b::Basis{p,Names},m::Monomial{p,Names};tag=b.notag,page=LASTPAGE) where {p,Names}
    if isa(tag,Monomial)
        tag = tag=>zero(GF{p})
    end
    if isempty(b.data) || b.data[end].v<m
        push!(b.data,(v=m,tag=tag,page=page))
        push!(b.lookup,m=>length(b.data))
    else
        pos = searchsortedfirst(b.data,(m,),by=first)
        insert!(b.data,pos,(v=m,tag=tag,page=page))
        push!(b.lookup,m=>pos)
        for (n,i)=b.lookup
            if i≥pos b.lookup[n] += 1 end
        end
    end
    b
end

function Base.deleteat!(b::Basis, where)
    deleteat!(b.data,where)
    if !isempty(where) && !all(==(false),where)
        empty!(b.lookup)
        for i=1:length(b.data)
            push!(b.lookup,b.data[i].v=>i)
        end
    end
    b
end

function gettag(b::Basis{p,Names},m::Monomial{p,Names}) where {p,Names}
    i,x = b[m]
    if i==0 || is_alive(x)
        nothing
    else
        x.tag
    end
end

function settag!(b::Basis{p,Names},i::Int,tag::Pair{Monomial{p,Names},GF{p}}=b.notag,page::Int=-1) where {p,Names}
    b.data[i] = (v=b.data[i].v,tag=tag,page=(page == -1 ? b.data[i].page : page))
end

################################################################
# Algebra
struct Algebra{p,Names}
    rules::Matrix{Union{Nothing,Vector{Tuple{GF{p},Gen{p,Names},Gen{p,Names}}}}}
    diff::Vector{Vector{Tuple{GF{p},Gen{p,Names},Gen{p,Names}}}}
    degree::Vector{Degree}
    basis::Dict{Degree,Basis{p,Names}}
    total::Ref{Int} # total degree to which we computed
    page::Ref{Int} # page number along spectral sequence
    dimension::Ref{Int} # record which sphere dimension we're looking at
end

const STABLE_DIMENSION = typemax(Int)

function Base.show(io::IO, Λ::Algebra{p,Names}) where {p,Names}
    print(io, "E", superscript_string(Λ.page[])," page of Λ algebra for S",Λ.dimension[]==STABLE_DIMENSION ? "ˢᵗᵃᵇˡᵉ" : superscript_string(Λ.dimension[])," over 𝔽",subscript_string(p)," at total degree≤",Λ.total[])
end
Base.show(io::IO, ::MIME"text/plain", Λ::Algebra) = show(io, Λ)

function Base.copy(Λ::Algebra{p,Names}) where {p,Names}
    Algebra{p,Names}(Λ.rules,Λ.diff,Λ.degree,Dict(k=>copy(b) for (k,b)=Λ.basis),Ref(Λ.total[]),Ref(Λ.page[]),Ref(Λ.dimension[]))
end

################################################################
# AlgebraElem
struct AlgebraElem{p,Names}
    parent::Algebra{p,Names}
    w::SortedDict{Monomial{p,Names},GF{p}}
#    w::Dict{Monomial{p,Names},GF{p}}
end
AlgebraElem(parent::Algebra{p,Names},w::SortedDict) where {p,Names} = AlgebraElem{p,Names}(parent,w)
AlgebraElem(parent::Algebra{p,Names},w::Monomial{p,Names},c::GF{p} = one(GF{p})) where {p,Names} = AlgebraElem(parent,SortedDict(w=>c))

Base.copy(x::AlgebraElem) = typeof(x)(x.parent,copy(x.w))
Base.:(==)(x::AlgebraElem{p,Names},y::AlgebraElem{p,Names}) where {p,Names} = x.w == y.w
Base.hash(x::AlgebraElem,h::UInt64) = hash(x.w,h)
Base.iterate(x::AlgebraElem,state...) = iterate(x.w,state...)
Base.length(x::AlgebraElem) = length(x.w)

function Base.show(io::IO, x::AlgebraElem)
    first = true
    for (k,v)=x
        iszero(v) && continue
        first || print(io,"+")
        isone(v) || print(io,v,"⋅")
        print(io,k)
        first = false
    end
    first && print(io,"𝟘")
end
#Base.show(io::IO, ::MIME"text/plain", x::AlgebraElem) = show(io, x)

Base.iszero(x::AlgebraElem) = all(iszero,values(x.w))
Base.zero(Λ::Algebra) = AlgebraElem(Λ,SortedDict())
Base.one(Λ::Algebra{p,Names}) where {p,Names} = AlgebraElem(Λ,SortedDict(Monomial{p,Names}()=>one(GF{p})))
Base.zero(x::AlgebraElem) = zero(x.parent)
Base.one(x::AlgebraElem) = one(x.parent)

function degree(Λ::Algebra{p,Names},m::Monomial{p,Names}) where {p,Names}
    d = (μ=0,λ=0,top=0)
    for i=m
        d = d + Λ.degree[i]
    end
    d
end

degree(Λ::Algebra{p,Names},g::Gen{p,Names}) where {p,Names} = Λ.degree[g]

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
    degree(x.parent,first(keys(x.w)))
end

function dimension(x::AlgebraElem)
    isempty(x.w) && return 0
    maximum(dimension(x.parent,β) for β=keys(x.w))
end

Base.mergewith(f,a::SortedDict{K,V},b::SortedDict{K,V}...) where {K,V} = (a = copy(a); mergewith!(f,a,b...); a)
    
function Base.:+(x::AlgebraElem{p,Names},y::AlgebraElem{p,Names}) where {p,Names}
    @assert x.parent == y.parent
    AlgebraElem(x.parent,mergewith(+,x.w,y.w))
end
function Base.:-(x::AlgebraElem{p,Names},y::AlgebraElem{p,Names}) where {p,Names}
    x + (-y)
end
function Base.:-(x::AlgebraElem)
    AlgebraElem(x.parent,SortedDict(k=>-v for (k,v)=x))
end
Base.:*(g::GF{p},x::AlgebraElem{p,Names}) where {p,Names} = iszero(g) ? zero(x) : AlgebraElem(x.parent,SortedDict(k=>g*v for (k,v)=x))
Base.:*(g::Integer,x::AlgebraElem{p,Names}) where {p,Names} = GF{p}(g)*x
Base.:^(x::AlgebraElem,y::Integer) = Base.power_by_squaring(x,y)

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

SET = 1
ADD = 0

"""add_monomial!(x::AlgebraElem,u,c,from=1)

add `c*copy(u)` to `x`.
 `u` is manipulated, to avoid allocation, but is preserved on exit.
`u` is guaranteed, on entry, to be reduced at least at positions `1:2`, `2:3`, ...,`from-1:from`.
"""
function add_monomial!(x::AlgebraElem{p,Names},u::Monomial{p,Names},c::GF{p}=one(GF{p}),from=1) where {p,Names}
    iszero(c) && return x
    
    rules = x.parent.rules
    for k=from:length(u)-1
        is_admissible_pair(u[k],u[k+1]) && continue
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
        global SET; SET += 1
        dict[copy(u)] = c # ... or make a copy of it
    else
        global ADD; ADD += 1
        newc = dict[stdu]+c
        if iszero(newc)
            delete!(dict,stdu)
        else
            dict[stdu] = newc
        end
    end    
    x
end
add_monomial!(x::AlgebraElem{p,Names},pair::Pair{Monomial{p,Names},GF{p}}) where {p,Names} = add_monomial!(x,pair...)

function leading_monomial(x::AlgebraElem{p,Names}) where {p,Names}
    last(x.w) # possible for SortedDict
end

"""tagger(Λ, m)

Searches through the bases of `Λ` to find the element n tagging `m`,
by selecting the minimal prefix of `m` such that the corresponding
suffix is in a basis; then return the prefix followed by the tagger
of the suffix.
"""
function deep_tagger(Λ::Algebra{p,Names}, m::Monomial, deg::Degree) where {p,Names}
    sign = one(GF{p})
    for i=1:length(m)-1
        deg = deg - Λ.degree[m[i]]
        if is_lambda(m[i])
            sign = -sign
        end
        γ = gettag(Λ.basis[deg],m[i+1:end])
        if γ≠nothing
            @assert !iszero(γ.second)
            return m[1:i]*γ.first => sign*γ.second
        end
    end
end

function reduce(Λ::Algebra{p,Names},x::AlgebraElem{p,Names}) where {p,Names}
    z = zero(x)
    z
    x
end

function Base.:*(x::AlgebraElem{p,Names},y::AlgebraElem{p,Names}) where {p,Names}
    Λ = x.parent
    @assert Λ == y.parent
    result = zero(x)
    for (xk,xv)=x, (yk,yv)=y
        @timeit_debug to "add_monomial(*)" add_monomial!(result,xk*yk,xv*yv)
    end
    reduce(Λ,result)
end

function add_differential!(result::AlgebraElem{p,Names}, u::Monomial{p,Names}, sign::GF{p}) where {p,Names}
    diff = result.parent.diff
    u = Gen{p,Names}(1)*u # prepare some space to store differential
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
add_differential!(result::AlgebraElem{p,Names}, m::Monomial{p,Names}) where {p,Names} = add_differential!(result, m, one(GF{p}))

differential(Λ::Algebra{p,Names}, m::Monomial{p,Names}) where {p,Names} = add_differential!(zero(Λ),m)
function differential(x::AlgebraElem)
    result = zero(x)
    for lm=x
        add_differential!(result,lm...)
    end
    result
end

"""Complete monomial `m` to a cycle
"""
function cycle(Λ::Algebra{p,Names},m::Monomial{p,Names}) where {p,Names}
    result = zero(Λ)
    add_monomial!(result,m)
    
    δ = differential(Λ,m)
    if !iszero(δ)
        δdegree = degree(δ)
        range = Λ[δdegree]
    end
    while !iszero(δ)
        m,c = leading_monomial(δ)
        j,y = range[m]
        if j==0
            @timeit_debug to "tagger" γ = deep_tagger(Λ,m,δdegree)
        else
            is_alive(y) && error("$m is not the leading monomial of a cycle")
            γ = y.tag
        end
        add_monomial!(result,γ.first,(-c)*γ.second)
        @timeit_debug to "add_differential" add_differential!(δ,γ.first,(-c)*γ.second)
    end
    return result
end

function differential_matrix(Λ::Algebra{p,Names},d::Degree) where {p,Names}
    source = Λ[d]
    range = Λ[d + (μ=0,λ=1,top=-1)]

    mat = zeros(GF{p},length(source),length(range))
    for (i,x)=enumerate(source)
        dx = differential(Λ,x.v)
        for (m,c)=dx
            j = range[m][1]
            mat[i,j] = c
        end
    end
    mat
end
    
"""Find a basis of cycles"""
function cycles(Λ::Algebra{p,Names},d::Degree) where {p,Names}
    
end

"""Find a basis of boundaries"""
function boundaries(Λ::Algebra{p,Names},d::Degree) where {p,Names}
    
end

function tag_basis!(Λ::Algebra{p,Names}, source::Basis{p,Names}, range::Basis{p,Names}) where {p,Names}
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
function prebasis(Λ::Algebra{p,Names},degree::Degree) where {p,Names}
    @assert all(≥(0),degree)

    b = Basis{p,Names}(degree)
        
    if degree.μ==degree.λ==0
        if degree.top==0
            push!(b, Monomial{p,Names}())
        end
        return b
    end

    for g::Gen{p,Names}=1:NGEN
        Λ.dimension[]<0 && !is_lambda(g) && continue

        newdegree = degree - Λ.degree[g]
        any(<(0),newdegree) && continue
        basis = Λ.basis[newdegree]
        for x=basis
            β = x.v
            is_admissible_product(g,β) || continue
            if !is_alive(x) && Λ.dimension[]≠0 # don't prune if dimension=0
                γ = x.tag.first
                is_admissible_product(g,γ) && continue
            end
            
            γ = g*β
            push!(b,γ)
        end
    end
    @assert issorted(b)
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
                @timeit_debug to "prebasis" Λ.basis[degree] = prebasis(Λ,degree)
            end
        end

        # we can do this in parallel, since all μ's operate independently
        #@threads # speedup seems to be only 50%
        for μ=0:Λ.total[]-1
            for λ=0:Λ.total[]-μ-1
                top = Λ.total[]-λ-μ
                source = (μ=μ,λ=λ,top=top)
                range = (μ=μ,λ=λ+1,top=top-1)
                @timeit_debug to "tag_basis!" tag_basis!(Λ,Λ.basis[source],Λ.basis[range])
            end
        end
    end
end

"""truncate!(Λ,dimension)

Remove all basis entries that have dimension too large
"""
function truncate!(Λ::Algebra{p,Names},max_dimension::Int) where {p,Names}
    foreach(values(Λ.basis)) do b
        deleteat!(b,[dimension(Λ,x.v)>max_dimension for x=b])
        # now manually erase tags that point to nowhere. could be optimized.
        for (i,x)=enumerate(b)
            if !is_alive(x) && dimension(Λ,x.tag.first)>max_dimension
                settag!(b,i,b.notag,LASTPAGE)
            end
        end
    end
    Λ.dimension[] = max_dimension
    Λ
end
Base.truncate(Λ::Algebra{p,Names},max_dimension::Int) where {p,Names} = truncate!(copy(Λ),max_dimension)

"""homology!(Λ)

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
    @assert iseven(count) # we delete twice, at source and range
    Λ.page[] += 1
    count÷2
end
homology(Λ::Algebra) = (L = copy(Λ); homology!(L); L)

function mark_differentials!(source::Vector{Basis{p,Names}},range::Vector{Basis{p,Names}},page::Int) where {p,Names}
    M = length(source)
    @assert M==length(range)

    # for all vectors in all source bases, in increasing order:
    # either give them a tag to range basis vector, & back, or mark them as untagged.

    count = 0
    survivors = Set{Monomial{p,Names}}() # we know these should not get killed
    
    for i=1:M
        for x=source[i]
            is_alive(x) || continue
            x.v∈survivors && continue
            taggees = [y.v for y=range[i] if is_alive(y)]
            for j=i+1:M
                x.v∈source[j] || setdiff!(taggees,[y.v for y=range[j]])
            end
            #!!! also check wrt Hopf maps
            if isempty(taggees)
                push!(survivors,x.v)
            else
                w = first(taggees)
                for j=i:M
                    l,_ = source[j][x.v]
                    m,_ = range[j][w]
                    if l>0 && m>0
                        count += 1
                        settag!(source[j],l,w=>zero(GF{p}),page)
                        settag!(range[j],m,x.v=>one(GF{p}),page)
                    end
                end
            end
        end
    end

    count
end

"""mark_differentials!(Λs::Vector{Algebra},Λ₊,dμ,dλ)

Add tags for higher differentials in each `Λ[n]`, with tridegree `(dμ,dλ,-1)`, that are
compatible with the suspension maps `Λ[n] → Λ[n+1]`.

Returns the number of new differentials.
"""
function mark_differentials!(Λ::Vector{Algebra{p,Names}}, d::Degree) where {p,Names}
    total = Λ[1].total[]
    page = Λ[1].page[]
    @assert all(L->L.total[]==total,Λ)
    @assert all(L->L.page[]==page,Λ)
    @assert d.top == -1
    @assert d.λ ≥ 1
    @assert d.μ ≥ 0
    
    count = 0

    for hom=total+1-d.μ-d.λ:-1:1
        for top=hom:total+1-hom-d.μ-d.λ
            for μ=0:hom
                sourcedegree = (μ=μ,λ=hom-μ,top=top)
                rangedegree = sourcedegree + d
                count += mark_differentials!([L.basis[sourcedegree] for L=Λ],[L.basis[rangedegree] for L=Λ],page)
            end
        end
    end

    count
end
mark_differentials!(Λ::Vector{Algebra{p,Names}}, μ::Int, λ::Int) where {p,Names} = mark_differentials!(Λ,(μ=μ,λ=λ,top=-1))

function Base.getindex(Λ::Algebra,degree::Degree)
    cache_basis!(Λ,degree.μ+degree.λ+degree.top)
    Λ.basis[degree]
end

Base.getindex(Λ::Algebra,i::Int,j::Int,top::Int) = Λ[(μ=i,λ=j,top=top)]

function Base.getindex(Λ::Algebra,::Colon,j::Int,top::Int)
    vcat((Λ.basis[(μ=μ,λ=j,top=top)].data for μ=0:Λ.total[]-j-top)...)
end

function Base.getindex(Λ::Algebra,ij::Int,top::Int)
    cache_basis!(Λ,ij+top)
    vcat((Λ.basis[(μ=μ,λ=ij-μ,top=top)].data for μ=0:ij)...)
end

function Base.getindex(Λ::Algebra,top::Int)
    vcat((Λ.basis[(μ=i,λ=j,top=top)].data for i=0:top for j=0:top-i)...)
end

function Base.getindex(Λ::Algebra{p,Names},m::Monomial{p,Names}) where {p,Names}
    d = degree(Λ,m)
    (d,Λ[d][m]...)
end
