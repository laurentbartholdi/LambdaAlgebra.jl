# unused

module ΛAlgebra
using Oscar
export lambda_algebra, homogeneous_basis, ∂, ∂matrix

Oscar.rank(m::SMat) = size(m,2)-nullspace(m)[1]

LOWERCASE_DIGIT = Dict('0'=>'₀','1'=>'₁','2'=>'₂','3'=>'₃','4'=>'₄','5'=>'₅','6'=>'₆','7'=>'₇','8'=>'₈','9'=>'₉')
lowercase_string(n) = prod(map(c->LOWERCASE_DIGIT[c],string(n)))
UPPERCASE_DIGIT = Dict('0'=>'⁰','1'=>'¹','2'=>'²','3'=>'³','4'=>'⁴','5'=>'⁵','6'=>'⁶','7'=>'⁷','8'=>'⁸','9'=>'⁹')
uppercase_string(n) = prod(map(c->UPPERCASE_DIGIT[c],string(n)))

const SHIFT = 1000 # all terms in monomial starting with 1000 are v's: 1 = λ₁ while 1000 = v₀

mutable struct Monomial
    c::FqFieldElem
    v::Vector{Int16}
    tags::Monomial
    taggedby::Monomial
    Monomial(c,v) = (m = new(); m.c = c; m.v = v; m)
end

struct LambdaAlgebra
    K::FqField
    R::FreeAssociativeAlgebra{FqFieldElem}
    rules::Matrix{Union{Nothing,Vector{Tuple{FqFieldElem,Int,Int}}}}
    diff::Vector{Vector{Tuple{FqFieldElem,Int,Int}}}
    degrees::Vector{Tuple{Int,Int,Int}}
    nlambdas::Int
    nvs::Int
    vbasis::Dict{Tuple{Int,Int},Vector{Pair{Vector{Int},Int}}}
    λbasis::Dict{Tuple{Int,Int},Vector{Vector{Int}}}
end

function Base.show(io::IO, Λ::LambdaAlgebra)
    print(io, "Λ algebra over ",Λ.K," with generators v₀:v",lowercase_string(Λ.nvs-1),",λ₁:λ",lowercase_string(Λ.nlambdas))
end

Oscar.characteristic(Λ::LambdaAlgebra) = characteristic(Λ.K)
    
struct LambdaAlgebraElem
    parent::LambdaAlgebra
    w::FreeAssociativeAlgebraElem{FqFieldElem}
end

function Base.show(io::IO, x::LambdaAlgebraElem)
    print(io, x.w)
end

Base.iszero(x::LambdaAlgebraElem) = iszero(x.w)
Base.isone(x::LambdaAlgebraElem) = isone(x.w)

function Oscar.is_homogeneous(x::LambdaAlgebraElem)
    iszero(x) && return true    
    d = lambda_degree(x.parent.degrees,x.w.exps[1])
    for i=1:x.w.length
        lambda_degree(x.parent.degrees,x.w.exps[i])==d || return false
    end
    return true
end

degree₊(u,v) = (u[1]+v[1],u[2]+v[2],u[3]+v[3])
degree₋(u,v) = (u[1]-v[1],u[2]-v[2],u[3]-v[3])

function lambda_degree(degrees,w)
    d = (0,0,0)
    for i=w
        d = degree₊(d,degrees[i])
    end
    d
end

function Oscar.degree(x::LambdaAlgebraElem)
    iszero(x) && error("0 element has no degree")
    d = lambda_degree(x.parent.degrees,x.w.exps[1])
    for i=1:x.w.length
        lambda_degree(x.parent.degrees,x.w.exps[i])==d || error("element is not homogeneous")
    end
    return d
end

function Base.:+(x::LambdaAlgebraElem,y::LambdaAlgebraElem)
    @assert x.parent == y.parent
    LambdaAlgebraElem(x.parent,x.w+y.w)
end
function Base.:-(x::LambdaAlgebraElem,y::LambdaAlgebraElem)
    @assert x.parent == y.parent
    LambdaAlgebraElem(x.parent,x.w-y.w)
end
function Base.:-(x::LambdaAlgebraElem)
    LambdaAlgebraElem(x.parent,-x.w)
end

function check(x::LambdaAlgebraElem)
    for u=x.w.exps[1:x.w.length]
        for k=1:length(u)-1
            if x.parent.rules[u[k],u[k+1]]≠nothing
                error("word $u should have been reduced")
            end
        end
    end
    return true
end

# return basis of v-part, as a vector of (vector=>weight), where weight(v_i)=p^i. the v_i's are decreasing.
function __v_basis(Λ::LambdaAlgebra,length::Int,degree::Int)
    get!(Λ.vbasis,(length,degree)) do
        if length==0
            return degree==0 ? [Int[] => 0] : Pair{Vector{Int},Int}[]
        end
        basis = Pair{Vector{Int},Int}[]
        for i=0:Λ.nvs-1
            pi = characteristic(Λ)^i
            geni = i + Λ.nlambdas+1
            Λ.degrees[geni][3] > degree && break
            for (v,weight)=__v_basis(Λ,length-1,degree-Λ.degrees[geni][3])
                if isempty(v) || v[1]≤geni
                    push!(basis,vcat(geni,v) => weight+pi)
                end
            end
        end
        basis
    end
end

# return basis of λ-part, as vector of vector. the condition is λ_{i+1} ≤ p*λ_i - 1
function __λ_basis(Λ::LambdaAlgebra,length::Int,degree::Int)
    get!(Λ.λbasis,(length,degree)) do
        if length==0
            return degree==0 ? [Int[]] : Vector{Int}[]
        end
        basis = Vector{Int}[]
        for i=1:Λ.nlambdas
            for v=__λ_basis(Λ,length-1,degree-Λ.degrees[i][3])
                if isempty(v) || v[1]≤characteristic(Λ)*i-1
                    push!(basis,vcat(i,v))
                end
            end
        end
        basis
    end
end

function __homogeneous_basis(Λ::LambdaAlgebra,degree::Tuple{Int,Int,Int},dimension)
    @assert isodd(dimension)
    dimension ÷= 2
    basis = Vector{Int}[]

    for d=0:degree[3]
        # combine vs of length/degree (degree[1],d) with λs of length/degree (degree[2],degree[3]-d)
        for (v,weight)=__v_basis(Λ,degree[1],d), w=__λ_basis(Λ,degree[2],degree[3]-d)
            if isempty(w) || w[1] ≤ dimension+weight
                push!(basis,vcat(v,w))
            end
        end
    end
    return basis
end

function homogeneous_basis_word_dict(Λ::LambdaAlgebra,degree::Tuple{Int,Int,Int},dimension=999)
    Dict(w=>i for (i,w)=enumerate(__homogeneous_basis(Λ,degree,dimension)))
end
    
function homogeneous_basis(Λ::LambdaAlgebra,degree::Tuple{Int,Int,Int},dimension=999)
    [LambdaAlgebraElem(Λ,AbstractAlgebra.Generic.FreeAssociativeAlgebraElem{FqFieldElem}(Λ.R,w)) for w=__homogeneous_basis(Λ,degree,dimension)]
end

# dict contains a polynomial, in the form monomial=>coeff; add c*copy(u) to it. u is preserved on exit.
function add_monomial!(dict,rules,u,c)
    stdu = getkey(dict,u,nothing)
    if stdu≠nothing
        dict[stdu] += c
        return
    end
    for k=1:length(u)-1
        fixes = rules[u[k],u[k+1]]
        if fixes≠nothing
            uk,uknext = u[k], u[k+1]
            for s=1:length(fixes)
                u[k] = fixes[s][2]
                u[k+1] = fixes[s][3]
                newc = fixes[s][1]*c
                add_monomial!(dict,rules,u,newc)
            end
            u[k],u[k+1] = uk,uknext
            return
        end
    end
    u = copy(u)
    if haskey(dict,u) dict[u] += c else dict[u] = c end
end
                          
function Base.:*(x::LambdaAlgebraElem,y::LambdaAlgebraElem)
    Λ = x.parent
    @assert Λ == y.parent

    result = Dict{Vector{Int},FqFieldElem}()
    for i in 1:x.w.length, j in 1:y.w.length
        add_monomial!(result,x.parent.rules,vcat(x.w.exps[i],y.w.exps[j]),x.w.coeffs[i]*y.w.coeffs[j])
    end
    zcoeffs = collect(values(result))
    zexps = collect(keys(result))
    w = AbstractAlgebra.Generic.FreeAssociativeAlgebraElem{FqFieldElem}(Λ.R, zcoeffs, zexps, length(zcoeffs))
    combine_like_terms!(sort_terms!(w))
    LambdaAlgebraElem(x.parent,w)
end

function ∂(x::LambdaAlgebraElem)
    result = Dict{Vector{Int},FqFieldElem}()
    Λ = x.parent
    for i in 1:x.w.length
        u = copy(x.w.exps[i])
        pushfirst!(u,0) # prepare some space to store differential
        sign = x.w.coeffs[i]
        for s=2:length(u) # u already has an extra letter to store the differential of a letter
            ui = u[s]
            for (c,x,y)=Λ.diff[ui]
                u[s-1],u[s] = x,y
                add_monomial!(result,Λ.rules,u,c*sign)
            end
            u[s-1] = ui
            if ui≤Λ.nlambdas
                sign = -sign
            end
        end                                                  
    end
    zcoeffs = collect(values(result))
    zexps = collect(keys(result))
    w = AbstractAlgebra.Generic.FreeAssociativeAlgebraElem{FqFieldElem}(Λ.R, zcoeffs, zexps, length(zcoeffs))
    combine_like_terms!(sort_terms!(w))
    LambdaAlgebraElem(x.parent,w)
end

# return the matrix of ∂, as a map from homogeneous_basis(Λ,degree) to homogeneous_basis(Λ,degree+(0,1,0)).
# the ith row of the matrix is the image of the ith basis vector.
function ∂matrix(Λ::LambdaAlgebra,degree::Tuple{Int,Int,Int},dimension=999)
    src = homogeneous_basis(Λ,degree,dimension)
    dstlookup = homogeneous_basis_word_dict(Λ,degree₊(degree,(0,1,0)),dimension)
                                                            
    """m = zero_matrix(Λ.K,length(src),length(dstlookup))
    for i=1:length(src)
        x = ∂(src[i])
        for s=1:x.w.length
            j = dstlookup[x.w.exps[s]]
            m[i,j] = x.w.coeffs[s]
        end
                                                            end"""
    m = sparse_matrix(Λ.K,0,length(dstlookup))
    for i=1:length(src)
        x = ∂(src[i])
        push!(m,sparse_row(Λ.K,[dstlookup[x.w.exps[i]] for i=1:x.w.length],x.w.coeffs[1:x.w.length]))
    end
    m
end

function lambda_algebra(p=3,n=100)
    K = GF(p)
    lambdas = Symbol[]
    k = 1; while (2p-2)*k ≤ n push!(lambdas,Symbol(:λ,lowercase_string(k))); k += 1 end
    shift = length(lambdas)
    vs = Symbol[]
    k = 0; while 2*p^k-1 ≤ n push!(vs,Symbol(:v,lowercase_string(k))); k += 1 end
    R,gens = free_associative_algebra(K,vcat(lambdas,vs))

    rules = Matrix{Union{Nothing,Vector{Tuple{FqFieldElem,Int,Int}}}}(undef,shift+length(vs),shift+length(vs))
    fill!(rules,nothing)
                                                                                                
    for i=1:length(vs), j=i+1:length(vs)
        rules[shift+i,shift+j] = [(one(K),shift+j,shift+i)]
    end
    for i=1:shift, j=0:length(vs)-1
        if j==0 || i+p^(j-1) > n
            rules[i,shift+1+j] = [(one(K),shift+1+j,i)]
        else
            rules[i,shift+1+j] = [(one(K),shift+1+j,i),(one(K),shift+j,i+p^(j-1))]
        end
    end
    for i=1:shift, k=0:shift-p*i
        rule = []
        for j=0:k-1
            push!(rule,(one(K)*(-1)^(j+1)*binomial((p-1)*(k-j)-1,j),i+k-j,p*i+j))
        end
        rules[i,p*i+k] = rule
    end
    
    diff = Vector{Tuple{FqFieldElem,Int,Int}}[]
    for i=1:shift
        push!(diff,[(one(K)*(-1)^(j+1)*binomial((p-1)*(i-j)-1,j),i-j,j) for j=1:i-1])
    end
    push!(diff,[])
    for j=1:length(vs)-1
        push!(diff,[(one(K),shift+j,p^(j-1))])
    end
    degrees = Tuple{Int,Int,Int}[]
    for i=1:shift
        push!(degrees,(0,1,(2p-2)*i))
    end
    for j=0:length(vs)-1
        push!(degrees,(1,0,2*p^j-1))
    end
                                            
    Λ = LambdaAlgebra(K,R,rules,diff,degrees,shift,length(vs),Dict(),Dict())
    Λ, [LambdaAlgebraElem(Λ,w) for w=gens[1:shift]],[LambdaAlgebraElem(Λ,w) for w=gens[shift+1:end]]
end

end

"""
L,lambdas,vs = ΛAlgebra.lambda_algebra()

baspre = ΛAlgebra.homogeneous_basis(L,(3,1,63))
bassrc = ΛAlgebra.homogeneous_basis(L,(3,2,63))
basdst = ΛAlgebra.homogeneous_basis(L,(3,3,63))
srcmat = ΛAlgebra.∂matrix(L,(3,2,63)) |> matrix
dstmat = ΛAlgebra.∂matrix(L,(3,1,63)) |> matrix
ker = kernel(srcmat)
srcmat

imgs = Dict(b=>(bb = ΛAlgebra.degree₋(b,(0,1,0)); bb∈bas ? solve(kers[b],basmat[bb]) : zero_matrix(L.K,0,size(kers[b],1))) for b=keys(kers))

"""

nothing
