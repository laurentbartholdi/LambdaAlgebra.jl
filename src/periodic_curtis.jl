# Context-sensitive Curtis reduction for the odd-primary periodic algebra.
# See docs/unstable-periodic.md for the cancellation criterion and proof.

const PV{p}=Dict{Monomial{p,LAMBDAV},GF{p}}
mutable struct PeriodicOps{p}
    algebra::Algebra{p,LAMBDAV}
    products::Dict{Monomial{p,LAMBDAV},PV{p}}
    diffs::Dict{Monomial{p,LAMBDAV},PV{p}}
    limit::Int
end
PeriodicOps(A::Algebra{p,LAMBDAV},limit::Int) where p=PeriodicOps{p}(A,Dict(),Dict(),limit)
function pvadd!(out::PV{p},m,c) where p
    z=get(out,m,zero(GF{p}))+c
    if iszero(z)
        delete!(out,m)
    else
        out[m]=z
    end
end
function pvsub!(out::PV{p},v::PV{p},c) where p
    for (m,a) in v
        pvadd!(out,m,-c*a)
    end
end
function genmul(O::PeriodicOps{p},g,m) where p
    gm=g*m
    if is_admissible_product(g,m)
        return PV{p}(gm=>one(GF{p}))
    end
    haskey(O.products,gm) && return O.products[gm]
    result=PV{p}()
    # Normalize the suffix first; cache only non-admissible joins.
    for (a,x,y) in O.algebra.rules[g,m[1]]
        for (n,c) in genmul(O,y,m[2:end]), (r,b) in genmul(O,x,n)
            pvadd!(result,r,a*c*b)
        end
    end
    length(O.products) >= O.limit && empty!(O.products)
    O.limit > 0 && (O.products[gm]=result)
    result
end
function pvdiff(O::PeriodicOps{p},m) where p
    haskey(O.diffs,m) && return O.diffs[m]
    out=PV{p}()
    if !isempty(m)
        g=m[1];tail=m[2:end]
        for (a,x,y) in O.algebra.diff[g]
            for (n,c) in genmul(O,y,tail), (r,b) in genmul(O,x,n)
                pvadd!(out,r,a*c*b)
            end
        end
        sign=is_lambda(g) ? -one(GF{p}) : one(GF{p})
        for (n,c) in pvdiff(O,tail), (r,b) in genmul(O,g,n)
            pvadd!(out,r,sign*c*b)
        end
    end
    length(O.diffs) >= O.limit && empty!(O.diffs)
    O.limit > 0 && (O.diffs[m]=out)
    out
end
function pvleading(x)
    m=maximum(keys(x))
    m=>x[m]
end

mutable struct PeriodicCurtis{p} <: AbstractCurtisCache
    algebra::Algebra{p,LAMBDAV}
    tables::Dict{NTuple{3,Int},Vector{Basis{p,LAMBDAV}}} # (μ, λ+top, sphere bound)
    prune::Bool
    omitted::Int
    rows::Dict{Tuple{Int,Int,Int,Monomial{p,LAMBDAV}},PV{p}} # explicit normalized boundaries
    ops::PeriodicOps{p}
    inherited_rows::Dict{Tuple{Int,Int,Int,Monomial{p,LAMBDAV}},PV{p}}
    published::Set{Tuple{Int,Int}}
    cache_limit::Int
end
function PeriodicCurtis(A::Algebra{p,LAMBDAV};prune=true,cache_limit=50000) where p
    cache_limit >= 0 || throw(ArgumentError("cache_limit must be nonnegative"))
    PeriodicCurtis{p}(A,Dict(),prune,0,Dict(),PeriodicOps(A,cache_limit),Dict(),Set(),cache_limit)
end
function copy_curtis_cache(E::PeriodicCurtis{p},A) where p
    PeriodicCurtis{p}(A,E.tables,E.prune,E.omitted,E.rows,E.ops,E.inherited_rows,copy(E.published),E.cache_limit)
end
function context_table(E::PeriodicCurtis{p},mu::Int,t::Int,D::Int) where p
    q=2p-2
    # Larger bounds admit exactly the same monomials in this diagonal.
    D=min(D,2div(t,q)+1)
    key=(mu,t,D)
    haskey(E.tables,key) && return E.tables[key]
    bs=Basis{p,LAMBDAV}[]
    for l=0:div(t,q)
        d=(μ=mu,λ=l,top=t-l)
        b=Basis{p,LAMBDAV}(d)
        if mu==l==t==0
            push!(b,Monomial{p,LAMBDAV}())
        else
            for gx in eachindex(E.algebra.degree)
                g=Gen{p,LAMBDAV}(gx)
                dg=E.algebra.degree[g]
                nd=d-dg
                any(<(0),nd) && continue
                if is_lambda(g)
                    mu==0 && 2index(g)+1<=D || continue
                    DD=2p*index(g)-1 # pure lambda suffix: next index ≤ p*i-1
                else
                    DD=D+dg.top+2 # removing v_i raises the suffix bound by 2p^i
                end
                nt=nd.λ+nd.top
                nt%q==0 && nd.λ<=div(nt,q) || continue
                cb=context_table(E,nd.μ,nt,DD)[nd.λ+1]
                for x in cb
                    is_admissible_product(g,x.v) || continue
                    if E.prune && !is_alive(x) && is_admissible_product(g,x.tag.first)
                        E.omitted+=1
                        continue
                    end
                    push!(b,g*x.v)
                end
            end
        end
        push!(bs,b)
        l>0 && context_tag!(E,bs[l],b,D)
    end
    E.tables[key]=bs
    bs
end
function context_deep(E::PeriodicCurtis{p},m::Monomial{p,LAMBDAV},deg::Degree,D::Int) where p
    sign=one(GF{p})
    for i=1:length(m)-1
        g=m[i]
        deg-=E.algebra.degree[g]
        if is_lambda(g)
            sign=-sign
            D=2p*index(g)-1
        else
            D+=E.algebra.degree[g].top+2
        end
        bs=context_table(E,deg.μ,deg.λ+deg.top,D)
        tag=gettag(bs[deg.λ+1],m[i+1:end])
        if tag!==nothing
            iszero(tag.second) && error("source suffix $m")
            return m[1:i]*tag.first=>sign*tag.second
        end
    end
    error("missing tag $m at $D $deg")
end
function context_row(E::PeriodicCurtis{p},m::Monomial{p,LAMBDAV},deg::Degree,D::Int) where p
    D=min(D,2div(deg.λ+deg.top,2p-2)+1)
    key=(deg.μ,deg.λ+deg.top,D,m)
    haskey(E.rows,key) && return E.rows[key]
    haskey(E.inherited_rows,key) && return E.inherited_rows[key]
    g=m[1]
    # With a cycle prefix there is no d(g)X correction.
    if isempty(E.algebra.diff[g])
        nd=deg-E.algebra.degree[g]
        DD=is_lambda(g) ? 2p*index(g)-1 : D+E.algebra.degree[g].top+2
        Y=context_row(E,m[2:end],nd,DD)
        return PV{p}(g*w=>a for (w,a) in Y)
    end
    # A stored tag records the leading source, not its completed preimage.
    # Eliminate larger pivots to reconstruct its normalized full boundary.
    tag=context_deep(E,m,deg,D)
    delta=PV{p}(m=>tag.second*c for (m,c) in pvdiff(E.ops,tag.first))
    while true
        n,c=pvleading(delta)
        if n==m
            isone(c) || error("bad lifted coefficient $m $c")
            break
        end
        isless(m,n) || error("missing leading target $m: $delta")
        Y=context_row(E,n,deg,D)
        pvsub!(delta,Y,c)
    end
    length(E.inherited_rows) >= E.cache_limit && empty!(E.inherited_rows)
    E.cache_limit > 0 && (E.inherited_rows[key]=delta)
    delta
end
function context_tag!(E::PeriodicCurtis{p},source,range,D) where p
    for (i,x) in enumerate(source)
        is_alive(x) || continue
        delta=copy(pvdiff(E.ops,x.v))
        while !isempty(delta)
            m,c=pvleading(delta)
            dimension(E.algebra,m)<=D || error("unstable differential $m from $(x.v) at $D")
            j,y=range[m]
            if j!=0 && is_alive(y)
                settag!(source,i,m=>zero(GF{p}),1)
                settag!(range,j,x.v=>inv(c),1)
                E.rows[(range.degree.μ,range.degree.λ+range.degree.top,D,m)]=PV{p}(m=>inv(c)*a for (m,a) in delta)
                break
            end
            Y=context_row(E,m,range.degree,D)
            pvsub!(delta,Y,c)
        end
    end
end

function cache_basis!(A::Algebra{p,LAMBDAV},total::Int,mu_limit::Int=-1) where p
    total < 0 && return A
    mu_limit >= -1 || throw(ArgumentError("μ_degree must be -1 or nonnegative"))
    total <= (2p-2)*(NGEN-NV-1) || throw(ArgumentError("requested total degree exceeds the available periodic generators"))
    if A.dimension == STABLE_DIMENSION
        invoke(cache_basis!,Tuple{Algebra,Int,Int},A,total,mu_limit)
        return A
    end
    A.page == 1 || throw(ArgumentError("extend the chain complex before taking homology"))
    A.curtis === nothing && (A.curtis=PeriodicCurtis(A))
    _publish_periodic!(A,A.curtis,total,mu_limit)
end
function _publish_periodic!(A,E::PeriodicCurtis{p},total,mu_limit) where p
    for t=0:2p-2:total, mu=0:(mu_limit < 0 ? total-t : min(mu_limit,total-t))
        (mu,t) in E.published && continue
        bs=context_table(E,mu,t,A.dimension)
        # The shared suffix tables must never be mutated by homology/truncate.
        for b in bs
            A.basis[b.degree]=copy(b)
        end
        push!(E.published,(mu,t))
    end
    A.total=max(A.total,total)
    A
end
function deep_tagger(A::Algebra{p,LAMBDAV},m::Monomial,deg::Degree) where p
    A.curtis === nothing && return invoke(deep_tagger,Tuple{Algebra,Monomial,Degree},A,m,deg)
    context_deep(A.curtis,m,deg,A.dimension)
end

function truncate!(A::Algebra{p,LAMBDAV},D::Int) where p
    isodd(D) && D>0 || throw(ArgumentError("the unstable periodic algebra requires a positive odd sphere dimension"))
    A.page == 1 || throw(ArgumentError("truncate the periodic chain complex before taking homology"))
    total=A.total
    mu_limit=maximum((d.μ for d in keys(A.basis));init=0)
    prune=A.curtis === nothing ? true : A.curtis.prune
    limit=A.curtis === nothing ? 50000 : A.curtis.cache_limit
    A.dimension=D
    empty!(A.basis)
    A.total=-1
    # Recompute the context: simply deleting stable tags loses unstable classes.
    A.curtis=PeriodicCurtis(A;prune,cache_limit=limit)
    cache_basis!(A,total,mu_limit)
end

"""curtis_stats(A)

Storage and cancellation counts for the context-sensitive periodic Curtis algorithm.
`stored` includes all auxiliary suffix contexts, not just the requested sphere.
"""
function curtis_stats(A::Algebra{p,LAMBDAV}) where p
    E=A.curtis
    E === nothing && return nothing
    (contexts=length(E.tables),
     stored=sum(length(b) for bs in values(E.tables) for b in bs;init=0),
     omitted=E.omitted,
     explicit_rows=length(E.rows),
     inherited_rows=length(E.inherited_rows),
     products=length(E.ops.products),
     differentials=length(E.ops.diffs))
end

function Base.truncate(A::Algebra{p,LAMBDAV},D::Int) where p
    L=copy(A)
    L.parent=nothing
    truncate!(L,D)
end
