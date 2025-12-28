# LambdaAlgebra

[![Build Status](https://github.com/laurentbartholdi/LambdaAlgebra.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/laurentbartholdi/LambdaAlgebra.jl/actions/workflows/CI.yml?query=branch%3Amain)

Computations in the lambda algebra, very primitive. Basic usage:

```julia
using LambdaAlgebra
L,ls,ms = lambda_algebra(p=3,dimension=0,top_degree=80)

Ls = [truncate(L,n) for n=1:50]
foreach(homology!,Ls) # the unstable lambda algebra, for S^n up to n=50

print_stem(LL,21:50,odd=true) # print the stem data

for k=1:15 # kill differentials by the "action principle"
    for λ=1:k
        μ = k-λ
        count = LambdaAlgebra.mark_differentials!(LL,μ,λ)
        @info "Removed $count arrows in degree ($μ,$λ,-1)"
    end
end

print_stem(LL,34,odd=true)

print_grid(LL[3],backend=:html) # print a page of the spectral sequence, in html

open("xxx.tex","w") do f; print_grid(f,LL[3],backend=:tikz) end # or in tex
```
