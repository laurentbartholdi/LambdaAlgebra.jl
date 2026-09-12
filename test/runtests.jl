using LambdaAlgebra
using Test
const LA = LambdaAlgebra

signature(A) = Dict(d=>count(LA.is_alive,b) for (d,b) in A.basis if any(LA.is_alive,b))
positive_signature(A) = filter(kv->kv.first.top>0,signature(A))

# Independent enumeration: choose the commuting v-prefix, then the lambda word
# from left to right. This does not call prebasis or any Curtis code.
function full_basis(A::LA.Algebra{p,LA.LAMBDAV},mu,t,l,D) where p
    M=LA.Monomial{p,LA.LAMBDAV}; G=LA.Gen{p,LA.LAMBDAV}
    out=M[]; word=G[]; q=2p-2
    function lambdas(left,weight,bound)
        if left==0
            weight==0 && push!(out,M(copy(word)))
            return
        end
        for i=1:min(bound,weight-left+1)
            push!(word,G(i))
            lambdas(left-1,weight-i,p*i-1)
            pop!(word)
        end
    end
    function polynomial(left,remaining,lo,weight)
        if left==0
            remaining%q==0 && lambdas(l,remaining÷q,(D-1)÷2+weight)
            return
        end
        for gx=lo:length(A.degree)
            g=G(gx)
            dt=A.degree[g].top
            dt<=remaining || break
            push!(word,g)
            polynomial(left-1,remaining-dt,gx,weight+p^LA.index(g))
            pop!(word)
        end
    end
    polynomial(mu,t,LA.NGEN-LA.NV,0)
    sort!(out)
end
function rank_mod_p(mat,p)
    a=copy(mat); r=0
    for j=1:size(a,2)
        k=findfirst(i->a[i,j]!=0,r+1:size(a,1))
        k===nothing && continue
        k+=r
        r+=1
        a[r,:],a[k,:]=a[k,:],a[r,:]
        a[r,:]=mod.(a[r,:]*invmod(a[r,j],p),p)
        for i=r+1:size(a,1)
            a[i,j]==0 && continue
            a[i,:]=mod.(a[i,:]-a[i,j]*a[r,:],p)
        end
        r==size(a,1) && break
    end
    r
end
function matrix_signature(A::LA.Algebra{p,LA.LAMBDAV},T,D) where p
    h=Dict{LA.Degree,Int}(); ops=LA.PeriodicOps(A,1000)
    for t=0:2p-2:T, mu=0:T-t
        bs=[full_basis(A,mu,t,l,D) for l=0:t÷(2p-2)]
        ranks=zeros(Int,length(bs)+1)
        for j=1:length(bs)-1
            lookup=Dict(m=>i for (i,m) in enumerate(bs[j+1]))
            mat=zeros(Int,length(bs[j]),length(bs[j+1]))
            for (i,m) in enumerate(bs[j])
                delta=differential(A,m)
                @test Dict(w=>c for (w,c) in delta if !iszero(c)) == LA.pvdiff(ops,m)
                @test iszero(differential(delta))
                for (w,c) in delta
                    iszero(c) && continue
                    @test haskey(lookup,w)
                    mat[i,lookup[w]]=Int(c)
                end
            end
            ranks[j+1]=rank_mod_p(mat,p)
        end
        for (j,b) in enumerate(bs)
            n=length(b)-ranks[j]-ranks[j+1]
            @test n>=0
            n>0 && (h[(μ=mu,λ=j-1,top=t-j+1)]=n)
        end
    end
    h
end

@testset "Unstable periodic Curtis reduction" begin
    @testset "Independent full matrices, p=$p, S^$D" for (p,D,T) in ((3,1,24),(3,3,28),(3,5,24),(5,3,32),(5,7,32))
        A,_,_=periodic_algebra(p=p,dimension=D,top_degree=T)
        @test signature(A)==matrix_signature(A,T,D)
        for (d,b) in A.basis, x in b
            LA.is_alive(x) || continue
            z=cycle(A,x.v)
            @test iszero(differential(z))
            @test dimension(z)<=D
            @test LA.leading_monomial(z).first==x.v
        end
    end
    @testset "Comparison with usual lambda-mu, p=$p, S^$D" for (p,D,T) in ((3,1,40),(3,3,64),(3,5,40),(5,3,40),(5,7,40))
        A,_,_=periodic_algebra(p=p,dimension=D,top_degree=T)
        B,_,_=lambda_algebra(p=p,dimension=D,top_degree=T)
        @test positive_signature(A)==positive_signature(B)
        # The current usual complex excludes words ending in mu; consequently
        # it omits the degree-zero tower present in the periodic complex.
        @test all(get(signature(A),(μ=k,λ=0,top=0),0)==1 for k=0:T)
    end
    @testset "Pruning and incremental computation" begin
        A,_,_=periodic_algebra(p=3,dimension=3,top_degree=48)
        B,_,_=periodic_algebra(p=3,dimension=3,top_degree=48,curtis=false)
        @test signature(A)==signature(B)
        @test LA.curtis_stats(A).omitted>0
        @test sum(length,values(A.basis)) < sum(length,values(B.basis))÷4
        C,_,_=periodic_algebra(p=3,dimension=3,top_degree=32,μ_degree=0)
        LA.cache_basis!(C,24,4) # increasing mu even though total decreases
        @test all(get(signature(C),d,0)==n for (d,n) in signature(A) if sum(d)<=24 && d.μ<=4)
        LA.cache_basis!(C,48)
        @test signature(C)==signature(A)
        before=signature(A)
        H=homology(A)
        @test signature(A)==before
        @test signature(H)==before
        @test H.parent===A
        @test all(LA.is_alive(x) for b in values(H.basis) for x in b)
        # Disabling disposable caches must change neither tags nor homology.
        N,_,_=periodic_algebra(p=3,dimension=3,top_degree=40,cache_limit=0)
        @test signature(N)==filter(kv->sum(kv.first)<=40,before)
    end
    @testset "Changing the sphere recomputes cancellations" begin
        A,_,_=periodic_algebra(p=3,top_degree=28)
        B=truncate(A,3)
        C,_,_=periodic_algebra(p=3,dimension=3,top_degree=28)
        @test signature(B)==signature(C)
        @test A.dimension==LA.STABLE_DIMENSION
        @test homology(B).parent===B
        @test_throws ArgumentError truncate(homology(B),1)
        truncate!(B,1)
        C1,_,_=periodic_algebra(p=3,dimension=1,top_degree=28)
        @test signature(B)==signature(C1)
    end
    @testset "A new smaller source after prefixing v0, pruning=$prune" for prune in (true,false)
        A,l,v=periodic_algebra(p=3,dimension=3,top_degree=28,μ_degree=4,curtis=prune)
        monomial(w)=LA.leading_monomial(w).first
        x=v[2]^3*l[2]*l[1]
        a=v[1]^3*l[5]*l[1]
        z=v[1]^3*l[3]*l[2]*l[1]
        vx=v[1]*x; va=v[1]*a; vz=v[1]*z
        @test dimension(a)==5
        @test dimension(va)==3
        # On S^3, x owns z. This tag must NOT be multiplied by v0.
        @test A[degree(x)][monomial(x)][2].tag.first==monomial(z)
        @test A[degree(z)][monomial(z)][2].tag==(monomial(x)=>LA.GF{3}(1))
        # Removing v0 changes the suffix context to S^5, where the smaller
        # source a is allowed and x survives instead of owning z.
        child=LA.context_table(A.curtis,3,24,5)
        @test LA.is_alive(child[3][monomial(x)][2])
        expected_tag=monomial(a)=>LA.GF{3}(-1)
        if prune
            @test child[4][monomial(z)][1]==0
            @test LA.context_deep(A.curtis,monomial(z),degree(z),5)==expected_tag
            @test A[degree(va)][monomial(va)][1]==0
            @test A[degree(vz)][monomial(vz)][1]==0
            @test LA.deep_tagger(A,monomial(vz),degree(vz))==(monomial(va)=>LA.GF{3}(-1))
        else
            @test child[4][monomial(z)][2].tag==expected_tag
            @test A[degree(vz)][monomial(vz)][2].tag==(monomial(va)=>LA.GF{3}(-1))
        end
        @test LA.is_alive(A[degree(vx)][monomial(vx)][2])
        expected_cycle=v[1]*v[2]^3*(l[2]*l[1]+l[1]*l[2]) +
            v[1]^4*(l[5]*l[1]+l[4]*l[2]-l[2]*l[4]) -
            v[1]^3*v[2]*l[2]*l[3]
        c=cycle(A,monomial(vx))
        @test c==expected_cycle
        @test iszero(differential(c))
        @test dimension(c)==3
        @test monomial(c)==monomial(vx)
        @test degree(c)==(μ=4,λ=2,top=22)
    end
    @testset "Generator capacity and input bounds" begin
        A,l,v=periodic_algebra(p=3,dimension=3)
        @test degree(l[125]).top==499
        @test length(v)==6
        @test_throws ArgumentError periodic_algebra(dimension=2)
        @test_throws ArgumentError periodic_algebra(dimension=0)
        @test_throws ArgumentError periodic_algebra(dimension=3,cache_limit=-1)
        @test_throws ArgumentError LA.cache_basis!(A,(2*3-2)*(LA.NGEN-LA.NV-1)+1)
    end
    @testset "Basis insertion indices" begin
        M=LA.Monomial{3,LA.LAMBDAV}; G=LA.Gen{3,LA.LAMBDAV}
        b=LA.Basis{3,LA.LAMBDAV}((μ=0,λ=1,top=0))
        for i in (3,1,2)
            push!(b,M([G(i)]))
        end
        @test all(b[x.v][1]==i for (i,x) in enumerate(b))
    end
end
