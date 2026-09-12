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
function pvsub!(out::PV{p},v,c) where p
    for (m,a) in v
        pvadd!(out,m,-c*a)
    end
end
function genmul(O::PeriodicOps{p},g,m) where p
    if is_admissible_product(g,m)
        return PV{p}(g*m=>one(GF{p}))
    elseif !is_lambda(g)
        # A v-generator only needs insertion in the sorted polynomial prefix.
        # No Adem expansion or intermediate dictionaries are involved.
        word=copy(m.v)
        pos=1
        while pos<=length(word) && !is_lambda(word[pos]) && word[pos]<g
            pos+=1
        end
        insert!(word,pos,g)
        return PV{p}(Monomial(word)=>one(GF{p}))
    end
    gm=g*m
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

# Even the smallest v-prefix has weight mu. Above this bound every lambda
# word in the diagonal is allowed, so all larger contexts are identical.
periodic_bound(mu,t,D,p)=min(D,max(1,2div(t,2p-2)+1-2mu))

mutable struct PeriodicCurtis{p} <: AbstractCurtisCache
    algebra::Algebra{p,LAMBDAV}
    tables::Dict{NTuple{3,Int},Vector{Basis{p,LAMBDAV}}} # (μ, λ+top, sphere bound)
    prune::Bool
    omitted::Int
    rows::Dict{Tuple{Int,Int,Int,WordID},WordPoly{p}}
    ops::PeriodicWords{p}
    inherited_rows::PeriodicCache{Tuple{Int,Int,Int,WordID},WordPoly{p}}
    lookups::Dict{NTuple{4,Int},Dict{WordID,Int}}
    published::Set{Tuple{Int,Int}}
    cache_limit::Int
end
function PeriodicCurtis(A::Algebra{p,LAMBDAV};prune=true,cache_limit=50000) where p
    cache_limit >= 0 || throw(ArgumentError("cache_limit must be nonnegative"))
    PeriodicCurtis{p}(A,Dict(),prune,0,Dict(),PeriodicWords(A,cache_limit),
        PeriodicCache{Tuple{Int,Int,Int,WordID},WordPoly{p}}(cache_limit),Dict(),Set(),cache_limit)
end
function copy_curtis_cache(E::PeriodicCurtis{p},A) where p
    PeriodicCurtis{p}(A,E.tables,E.prune,E.omitted,E.rows,E.ops,E.inherited_rows,E.lookups,copy(E.published),E.cache_limit)
end
function context_table(E::PeriodicCurtis{p},mu::Int,t::Int,D::Int) where p
    q=2p-2
    # Larger bounds admit exactly the same monomials in this diagonal.
    D=periodic_bound(mu,t,D,p)
    key=(mu,t,D)
    haskey(E.tables,key) && return E.tables[key]
    bs=Basis{p,LAMBDAV}[]
    for l=0:div(t,q)
        d=(μ=mu,λ=l,top=t-l)
        b=Basis{p,LAMBDAV}(d)
        if mu==l==t==0
            push!(b,Monomial{p,LAMBDAV}())
        else
            generators=mu==0 ? (1:min(div(t,q),(D-1)÷2)) : ((NGEN-NV):length(E.algebra.degree))
            for gx in generators
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
        E.lookups[(mu,t,D,l)]=Dict(wordid(E.ops,x.v)=>i for (i,x) in enumerate(b))
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
function worddeep(E::PeriodicCurtis{p},m::WordID,deg::Degree,D::Int) where p
    O=E.ops; prefix=Int16[]; sign=one(GF{p})
    while O.tails[m]!=0
        g=O.heads[m]; m=wordtail(O,m)
        m==1 && break
        push!(prefix,g)
        dg=E.algebra.degree[g]; deg-=dg
        if g<NGEN-NV
            sign=-sign; D=2p*g-1
        else
            D+=dg.top+2
        end
        t=deg.λ+deg.top; D=periodic_bound(deg.μ,t,D,p)
        bs=context_table(E,deg.μ,t,D)
        j=get(E.lookups[(deg.μ,t,D,deg.λ)],m,0)
        if j!=0
            x=bs[deg.λ+1][j]
            if !is_alive(x)
                iszero(x.tag.second) && error("source suffix")
                source=wordid(O,x.tag.first)
                for h in Iterators.reverse(prefix)
                    source=wordcons(O,h,source)
                end
                return source=>sign*x.tag.second
            end
        end
    end
    error("missing implied word tag")
end
function context_wordrow(E::PeriodicCurtis{p},m::WordID,deg::Degree,D::Int) where p
    D=periodic_bound(deg.μ,deg.λ+deg.top,D,p)
    key=(deg.μ,deg.λ+deg.top,D,m)
    haskey(E.rows,key) && return E.rows[key]
    cached=get(E.inherited_rows,key,nothing)
    cached!==nothing && return cached
    O=E.ops; g=O.heads[m]
    # With a cycle prefix there is no d(g)X correction.
    if isempty(E.algebra.diff[g])
        nd=deg-E.algebra.degree[g]
        DD=g<NGEN-NV ? 2p*g-1 : D+E.algebra.degree[g].top+2
        Y=context_wordrow(E,wordtail(O,m),nd,DD)
        return wordprefix(O,g,Y)
    end
    # A stored tag records the leading source, not its completed preimage.
    # Eliminate larger pivots to reconstruct its normalized full boundary.
    tag=worddeep(E,m,deg,D)
    delta=WordPoly{p}(m=>tag.second*c for (m,c) in worddiff(O,tag.first))
    while true
        n,c=wordleading(O,delta)
        if n==m
            isone(c) || error("bad lifted coefficient $m $c")
            break
        end
        wordless(O,m,n) || error("missing leading target $m: $delta")
        Y=context_wordrow(E,n,deg,D)
        wordsub!(delta,Y,c)
    end
    E.inherited_rows[key]=delta
    delta
end
function context_row(E::PeriodicCurtis{p},m::Monomial{p,LAMBDAV},deg::Degree,D::Int) where p
    PV{p}(publicword(E.ops,w)=>c for (w,c) in context_wordrow(E,wordid(E.ops,m),deg,D))
end
function context_tag!(E::PeriodicCurtis{p},source,range,D) where p
    O=E.ops
    lookup=E.lookups[(range.degree.μ,range.degree.λ+range.degree.top,D,range.degree.λ)]
    for (i,x) in enumerate(source)
        is_alive(x) || continue
        delta=copy(worddiff(O,wordid(O,x.v)))
        while !isempty(delta)
            m,c=wordleading(O,delta)
            worddimension(O,m)<=D || error("unstable differential $m from $(x.v) at $D")
            j=get(lookup,m,0); y=j==0 ? nothing : range[j]
            if j!=0 && is_alive(y)
                settag!(source,i,publicword(O,m)=>zero(GF{p}),1)
                settag!(range,j,x.v=>inv(c),1)
                E.rows[(range.degree.μ,range.degree.λ+range.degree.top,D,m)]=WordPoly{p}(m=>inv(c)*a for (m,a) in delta)
                break
            end
            Y=context_wordrow(E,m,range.degree,D)
            wordsub!(delta,Y,c)
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
        length(E.ops.nodes)>=E.ops.next_collection && collect_periodic_words!(E)
    end
    A.total=max(A.total,total)
    A
end

# Run only between completed root diagonals: all live internal IDs are then
# reachable from these tables and caches. Freed IDs retain no external API.
function collect_periodic_words!(E::PeriodicCurtis)
    O=E.ops; marked=falses(length(O.heads)); marked[1]=true
    function mark(id)
        while id!=0 && !marked[id]
            marked[id]=true
            id=O.tails[id]
        end
    end
    for id in keys(O.public_words)
        mark(id)
    end
    for lookup in values(E.lookups), id in keys(lookup)
        mark(id)
    end
    for (key,row) in E.rows
        mark(key[4])
        foreach(mark,keys(row))
    end
    for cache in (E.inherited_rows.recent,E.inherited_rows.previous), (key,row) in cache
        mark(key[4])
        foreach(mark,keys(row))
    end
    for cache in (O.products.recent,O.products.previous,O.diffs.recent,O.diffs.previous), (id,row) in cache
        mark(id)
        foreach(mark,keys(row))
    end
    freed=0
    for id=2:length(O.heads)
        if O.heads[id]!=0 && !marked[id]
            delete!(O.nodes,(O.heads[id],O.runs[id],O.tails[id]))
            O.heads[id]=0; O.runs[id]=0; O.tails[id]=0
            push!(O.free,WordID(id))
            freed+=1
        end
    end
    O.next_collection=max(200000,2length(O.nodes))
    freed
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
     words=length(E.ops.nodes)+1,
     word_slots=length(E.ops.heads),
     products=length(E.ops.products),
     differentials=length(E.ops.diffs))
end

function Base.truncate(A::Algebra{p,LAMBDAV},D::Int) where p
    L=copy(A)
    L.parent=nothing
    truncate!(L,D)
end
