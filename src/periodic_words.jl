# Interned, run-length encoded words for the private periodic reducer.
# A node is (first generator, multiplicity, remaining word); ID 1 is empty.
const WordID=UInt32
const WordPoly{p}=Dict{WordID,GF{p}}
# Two generations retain recently reused entries when capacity is reached.
# Values are immutable to callers; eviction only causes recomputation.
mutable struct PeriodicCache{K,V}
    recent::Dict{K,V}
    previous::Dict{K,V}
    capacity::Int
end
PeriodicCache{K,V}(limit) where {K,V}=PeriodicCache{K,V}(Dict(),Dict(),limit)
Base.length(c::PeriodicCache)=length(c.recent)+length(c.previous)
function Base.setindex!(c::PeriodicCache,v,k)
    c.capacity==0 && return v
    if !haskey(c.recent,k) && length(c.recent)>=max(1,c.capacity÷2)
        c.previous=c.capacity==1 ? empty(c.recent) : c.recent
        c.recent=empty(c.recent)
    end
    delete!(c.previous,k)
    c.recent[k]=v
end
function Base.get(c::PeriodicCache,k,default)
    v=get(c.recent,k,nothing)
    v!==nothing && return v
    v=get(c.previous,k,nothing)
    v===nothing && return default
    delete!(c.previous,k)
    c[k]=v
    v
end
mutable struct PeriodicWords{p}
    algebra::Algebra{p,LAMBDAV}
    heads::Vector{Int16}
    runs::Vector{Int16}
    tails::Vector{WordID}
    nodes::Dict{Tuple{Int16,Int16,WordID},WordID}
    free::Vector{WordID}
    next_collection::Int
    public_words::Dict{WordID,Monomial{p,LAMBDAV}}
    products::PeriodicCache{WordID,WordPoly{p}}
    diffs::PeriodicCache{WordID,WordPoly{p}}
    binomials::Vector{Vector{GF{p}}}
end
function PeriodicWords(A::Algebra{p,LAMBDAV},limit) where p
    PeriodicWords{p}(A,[0],[0],[0],Dict(),[],200000,Dict(WordID(1)=>Monomial{p,LAMBDAV}()),
        PeriodicCache{WordID,WordPoly{p}}(limit),PeriodicCache{WordID,WordPoly{p}}(limit),[[one(GF{p})]])
end
function wordbinomials(O::PeriodicWords{p},n) where p
    while length(O.binomials)<=n
        previous=O.binomials[end]; row=ones(GF{p},length(previous)+1)
        for k=2:length(previous)
            row[k]=GF{p}(Int(previous[k-1])+Int(previous[k]))
        end
        push!(O.binomials,row)
    end
    O.binomials[n+1]
end
function wordnode(O::PeriodicWords,g::Int16,n::Int,tail::WordID)
    n==0 && return tail
    if O.heads[tail]==g
        n+=O.runs[tail]
        tail=O.tails[tail]
    end
    key=(g,Int16(n),tail)
    id=get(O.nodes,key,WordID(0))
    id!=0 && return id
    if isempty(O.free)
        id=WordID(length(O.heads)+1)
        push!(O.heads,g); push!(O.runs,Int16(n)); push!(O.tails,tail)
    else
        id=pop!(O.free)
        O.heads[id]=g; O.runs[id]=Int16(n); O.tails[id]=tail
    end
    O.nodes[key]=id
    id
end
wordcons(O,g::Int16,tail::WordID)=wordnode(O,g,1,tail)
function wordtail(O,id::WordID)
    O.runs[id]==1 ? O.tails[id] : wordnode(O,O.heads[id],Int(O.runs[id])-1,O.tails[id])
end
function wordcat(O,a::WordID,b::WordID)
    a==1 && return b
    wordnode(O,O.heads[a],Int(O.runs[a]),wordcat(O,O.tails[a],b))
end
function wordid(O::PeriodicWords{p},m::Monomial{p,LAMBDAV}) where p
    id=WordID(1)
    for g in Iterators.reverse(m.v)
        id=wordcons(O,g.x,id)
    end
    O.public_words[id]=m
    id
end
function publicword(O::PeriodicWords{p},id::WordID) where p
    haskey(O.public_words,id) && return O.public_words[id]
    word=Gen{p,LAMBDAV}[]
    node=id
    while node!=1
        g=Gen{p,LAMBDAV}(O.heads[node])
        for _=1:O.runs[node]
            push!(word,g)
        end
        node=O.tails[node]
    end
    m=Monomial(word)
    O.public_words[id]=m
    m
end
function wordless(O,a::WordID,b::WordID)
    while a!=b
        ga=O.heads[a]; gb=O.heads[b]
        ga!=gb && return ga<gb
        na=O.runs[a]; nb=O.runs[b]
        if na!=nb
            return na<nb ? O.heads[O.tails[a]]<ga : ga<O.heads[O.tails[b]]
        end
        a=O.tails[a]; b=O.tails[b]
    end
    false
end
function wordleading(O,v)
    m=WordID(1)
    for n in keys(v)
        wordless(O,m,n) && (m=n)
    end
    m=>v[m]
end
function worddimension(O,m::WordID)
    shift=0
    while m!=1
        g=O.heads[m]
        g<NGEN-NV && return max(0,2g+1-shift)
        shift+=Int(O.runs[m])*(O.algebra.degree[g].top+2)
        m=O.tails[m]
    end
    0
end
function wordadd!(out::WordPoly{p},m::WordID,c) where p
    z=get(out,m,zero(GF{p}))+c
    iszero(z) ? delete!(out,m) : (out[m]=z)
    out
end
function wordsub!(out::WordPoly,v,c)
    for (m,a) in v
        wordadd!(out,m,-c*a)
    end
    out
end
function wordv(O,g,m::WordID)
    h=O.heads[m]
    (h<NGEN-NV || h>=g) && return wordcons(O,g,m)
    wordnode(O,h,Int(O.runs[m]),wordv(O,g,O.tails[m]))
end
function word_admissible(O::PeriodicWords{p},g,h) where p
    h==0 || (g<NGEN-NV ? h<NGEN-NV && h<=p*g-1 : h<NGEN-NV || h>=g)
end
function wordmuladd!(out::WordPoly{p},O::PeriodicWords{p},g::Int16,m::WordID,c) where p
    if word_admissible(O,g,O.heads[m])
        wordadd!(out,wordcons(O,g,m),c)
    elseif g>=NGEN-NV
        wordadd!(out,wordv(O,g,m),c)
    else
        for (n,a) in wordproduct(O,g,m)
            wordadd!(out,n,c*a)
        end
    end
    out
end
function wordproduct(O::PeriodicWords{p},g::Int16,m::WordID) where p
    gm=wordcons(O,g,m)
    word_admissible(O,g,O.heads[m]) && return WordPoly{p}(gm=>one(GF{p}))
    g>=NGEN-NV && return WordPoly{p}(wordv(O,g,m)=>one(GF{p}))
    cached=get(O.products,gm,nothing)
    cached!==nothing && return cached
    result=WordPoly{p}(); h=O.heads[m]
    if h>=NGEN-NV
        # lambda_a v_j^n = sum binomial(n,k) v_{j-1}^k v_j^{n-k}
        #                         lambda_{a+k*p^(j-1)} (j>0).
        n=Int(O.runs[m]); tail=O.tails[m]
        j=Int(h)-(NGEN-NV)
        if j==0
            part=WordPoly{p}()
            wordmuladd!(part,O,g,tail,one(GF{p}))
            for (w,c) in part
                wordadd!(result,wordnode(O,h,n,w),c)
            end
        else
            coefficients=wordbinomials(O,n); shift=p^(j-1)
            for k=0:n
                a=coefficients[k+1]; iszero(a) && continue
                idx=g+k*shift
                idx<NGEN-NV || error("periodic product exceeds generator capacity")
                part=WordPoly{p}()
                wordmuladd!(part,O,Int16(idx),tail,a)
                for (w,c) in part
                    w=wordnode(O,h,n-k,w)
                    w=wordnode(O,Int16(h-1),k,w)
                    wordadd!(result,w,c)
                end
            end
        end
    else
        tail=wordtail(O,m)
        for (a,x,y) in O.algebra.rules[g,h]
            if word_admissible(O,y.x,O.heads[tail])
                wordmuladd!(result,O,x.x,wordcons(O,y.x,tail),a)
            else
                for (n,c) in wordproduct(O,y.x,tail)
                    wordmuladd!(result,O,x.x,n,a*c)
                end
            end
        end
    end
    O.products[gm]=result
    result
end
function worddiff(O::PeriodicWords{p},m::WordID) where p
    cached=get(O.diffs,m,nothing)
    cached!==nothing && return cached
    out=WordPoly{p}()
    if m!=1
        g=O.heads[m]; dg=O.algebra.diff[g]
        if isempty(dg)
            # Differentiate beyond an entire run of cycles in one step.
            tail=O.tails[m]; n=Int(O.runs[m])
            sign=g<NGEN-NV && isodd(n) ? -one(GF{p}) : one(GF{p})
            for (w,c) in worddiff(O,tail)
                word_admissible(O,g,O.heads[w]) || error("cycle prefix lost admissibility")
                wordadd!(out,wordnode(O,g,n,w),sign*c)
            end
        elseif g>NGEN-NV
            # Differentiate the whole commuting power using binomial
            # coefficients in the ground field, before expanding any terms.
            n=Int(O.runs[m]); tail=O.tails[m]
            coefficients=wordbinomials(O,n)
            shift=p^(Int(g)-(NGEN-NV)-1)
            for k=1:n
                a=coefficients[k+1]; iszero(a) && continue
                idx=k*shift
                idx<NGEN-NV || error("periodic differential exceeds generator capacity")
                part=WordPoly{p}()
                wordmuladd!(part,O,Int16(idx),tail,a)
                for (w,c) in part
                    w=wordnode(O,g,n-k,w)
                    w=wordnode(O,Int16(g-1),k,w)
                    wordadd!(out,w,c)
                end
            end
            for (w,c) in worddiff(O,tail)
                word_admissible(O,g,O.heads[w]) || error("v power lost admissibility")
                wordadd!(out,wordnode(O,g,n,w),c)
            end
        else
            tail=wordtail(O,m)
            for (a,x,y) in dg
                if word_admissible(O,y.x,O.heads[tail])
                    wordmuladd!(out,O,x.x,wordcons(O,y.x,tail),a)
                else
                    for (w,c) in wordproduct(O,y.x,tail)
                        wordmuladd!(out,O,x.x,w,a*c)
                    end
                end
            end
            sign=g<NGEN-NV ? -one(GF{p}) : one(GF{p})
            for (w,c) in worddiff(O,tail)
                wordmuladd!(out,O,g,w,sign*c)
            end
        end
    end
    O.diffs[m]=out
    out
end
struct WordRow{p}
    words::PeriodicWords{p}
    prefix::WordID
    terms::WordPoly{p}
end
Base.length(row::WordRow)=length(row.terms)
Base.eltype(::Type{WordRow{p}}) where p=Pair{WordID,GF{p}}
function Base.iterate(row::WordRow,state...)
    item=iterate(row.terms,state...)
    item===nothing && return nothing
    (m,c),next=item
    (wordcat(row.words,row.prefix,m)=>c),next
end
wordprefix(O,g::Int16,Y::WordPoly)=WordRow(O,wordcons(O,g,WordID(1)),Y)
wordprefix(O,g::Int16,Y::WordRow)=WordRow(O,wordcons(O,g,Y.prefix),Y.terms)
