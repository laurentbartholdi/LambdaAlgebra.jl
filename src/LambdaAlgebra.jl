module LambdaAlgebra

using TimerOutputs, PrettyTables, Base.Threads, DataStructures

export periodic_algebra, lambda_algebra
export degree, is_homogeneous, dimension, differential, cycle
export homology, homology!, truncate!, mark_differentials!
export print_curtis_table, print_grid, print_stem

const to = TimerOutput()

SUBSCRIPT_DIGIT = Dict('0'=>'₀','1'=>'₁','2'=>'₂','3'=>'₃','4'=>'₄','5'=>'₅','6'=>'₆','7'=>'₇','8'=>'₈','9'=>'₉','-'=>'₋','+'=>'⁺')
subscript_string(n) = prod(map(c->SUBSCRIPT_DIGIT[c],string(n)))
SUPERSCRIPT_DIGIT = Dict('0'=>'⁰','1'=>'¹','2'=>'²','3'=>'³','4'=>'⁴','5'=>'⁵','6'=>'⁶','7'=>'⁷','8'=>'⁸','9'=>'⁹','-'=>'⁻','+'=>'⁺')
superscript_string(n) = prod(map(c->SUPERSCRIPT_DIGIT[c],string(n)))

include("common.jl")

include("print.jl")

################################################################
# the periodic lambda algebra, with generators vᵢ,λᵢ
include("periodic.jl")

################################################################
# the classical lambda algebra, with generators λᵢ,μᵢ
include("lambdamu.jl")

################################################################
# the classical lambda algebra, with generators λᵢ, for p=2
include("lambda.jl")

# Wikipedia, https://en.wikipedia.org/wiki/Homotopy_groups_of_spheres

S³(q) = 1 + 2q^4 + 2q^5 + 12q^6 + 2q^7 + 2q^8 + 3q^9 + 15q^10 + 2q^11 + (2*2)q^12 + (12*2)q^13 + (84*2*2)q^14 + (2*2)q^15 + 6q^16 + 30q^17 + 30q^18 + (6*2)q^19 + (12*2*2)q^20 + (12*2*2)q^21 + (132*2)q^22

# Toda pp. 73-74
S³₃modA(q) = q^13 + # b
    q^16 + # ab
    q^17 + # AB
    q^23 + # b²
    q^24 + # B²
    q^26 + # ab²
    q^27 + # AB²
    q^32 + # ab₂
    2*q^33 + # AB₂,b³
    q^34 + # B³
    q^36 + # ab³
    q^37 + # AB³
    q^39 + # bb₂
    q^40 + # BB₂
    q^42 + # abb₂
    2*q^43 + # ABB₂,b⁴
    q^44 + # B⁴
    2*q^49 + # AE₂,b²b₂
    2*q^50 + # B²B₂,be'
    q^51 + # BE'
    q^52 + # ab²b₂
    2*q^53 + # AB²B₂,b⁵
    q^54 + # B⁵
    q^55 + # b₂²
    q^56 + # B₂²
    q^58 + # ab₂²
    2*q^59 + # AB₂²,b³b₂
    2*q^60 + # B³B₂,U
    q^61 + # B²E'
    q^62 + # ab³b₂
    q^63 + # AB³B₂
    q^65 + # bb₂² 
    q^66 + # BB₂²
    q^68 + # abb₂²
    q^69 + # ABB₂²
    q^75 + # b²b₂²
    q^76 + # B²B₂²
    q^78 + # m
    2*q^79 # AB²B₂²,???

S³₃(q) = q^3 + (1+q)*q^6/(1-q^4) + S³modA(q)

#= Toda pp. 60-61
    q^3 + q^6 + q^7 + q^10 +
    q^11 + q^13 + q^14 + q^15 + q^16 + q^17 + q^18 + q^19 +
    q^22 + 2*q^23 + q^24 + 2*q^26 + 2*q^27 + q^30 +
    q^31 + q^32 + 2*q^33 + 2*q^34 + q^35 + q^36 + q^37 + q^38 + 2*q^39 + q^40 +
    2*q^42 + 3*q^43 + q^44 + q^46 + q^47 + 2*q^49 + 3*q^50 +
    2*q^51 + q^52 + 2*q^53 + 2*q^54 + 2*q^55
=#

#================================================================ experiments

using LLLplus
function findrels(v,ϵ = sqrt(eps(real(eltype(v)))))
    n = length(v)
    m = maximum(abs.(v))
    mat = zeros(eltype(v),n+1,n)
    for i=1:n
        mat[i,i] = m*ϵ
        mat[n+1,i] = v[i]
    end
    LLLplus.lll(mat)[2]
end

const BigComplex = Complex{BigFloat}

Δ(q::T) where T = (T(2)*pi)^12*q*prod((T(1)-q^n) for n=1:1000)^24

mat = setprecision(2000) do
    z = BigComplex(0.7+0.6im)
    F = S³
#    F = Δ
    m = [1 0;8 1]
    findrels([BigComplex(pi)*2im,log(m[2,1]*z+m[2,2]),log(F(exp(BigComplex(2im)*pi*z))),log(F(exp(BigComplex(2im)*pi*(m[1,1]*z+m[1,2])/(m[2,1]*z+m[2,2]))))],1e-15)
end

================================================================#

end

#= products

1   (μ = 0, λ = 0, top = 0)  => [𝟙]
l   (μ = 0, λ = 1, top = 3)  => [λ₁]
ll  (μ = 0, λ = 2, top = 6)  => [λ₁²]
m   (μ = 0, λ = 2, top = 10) => [λ₁λ₂+λ₂λ₁]
lm  (μ = 0, λ = 3, top = 13) => [λ₁(λ₁λ₂+λ₂λ₁)]
n   (μ = 0, λ = 3, top = 21) => [λ₁λ₂λ₃]
llm (μ = 0, λ = 4, top = 16) => [λ₁²(λ₁λ₂+λ₂λ₁)]
o   (μ = 0, λ = 4, top = 24) => [λ₁²λ₂λ₃]
llmm(μ = 0, λ = 5, top = 23) => [λ₁²(λ₁λ₂+λ₂λ₁)²]
q   (μ = 1, λ = 2, top = 10) => [λ₁μ₁λ₁]
r   (μ = 1, λ = 2, top = 14) => [2⋅λ₁μ₁λ₂+λ₁μ₂λ₁]
lq  (μ = 1, λ = 3, top = 13) => [λ₁²μ₁λ₁]
s   (μ = 1, λ = 3, top = 17) => [λ₁²μ₁λ₂+λ₁μ₁λ₁λ₂+λ₁μ₁λ₂λ₁+λ₁λ₂μ₁λ₁]
t   (μ = 1, λ = 3, top = 25) => [λ₁μ₁λ₂λ₃]
u   (μ = 1, λ = 4, top = 20) => [λ₁μ₁λ₁(λ₁λ₂+λ₂λ₁)]
v   (μ = 1, λ = 4, top = 24) => [(2⋅λ₁μ₁λ₂+λ₁μ₂λ₁)(λ₁λ₂+λ₂λ₁)]
lu  (μ = 1, λ = 5, top = 23) => [λ₁²μ₁λ₁(λ₁λ₂+λ₂λ₁)]
w   (μ = 2, λ = 2, top = 14) => [λ₁μ₁²λ₁]
    (μ = 2, λ = 3, top = 21) => [λ₁²μ₁²λ₂+2⋅λ₁²μ₁μ₂λ₁+λ₁²μ₂μ₁λ₁+λ₁μ₁λ₁μ₁λ₂+λ₁μ₁²λ₁λ₂+λ₁μ₁²λ₂λ₁+λ₁μ₂λ₁μ₁λ₁]
  (μ = 2, λ = 4, top = 24) => [λ₁μ₁²λ₁(λ₁λ₂+λ₂λ₁)]
  (μ = 3, λ = 2, top = 18) => [λ₁μ₁³λ₁]
  (μ = 4, λ = 2, top = 22) => [λ₁μ₁⁴λ₁]

=#

#= a class

function xx(n)
    μ₁,μ₂,μ₃ = ms[2:4]
    λ₁,λ₂,λ₃ = ls[1:3]
    μ₁^(3n)*λ₃ + sum(i*μ₁^(3n-1-i)*μ₂*μ₁^i*λ₂ for i=0:3n-1) + sum((1+i+j)*(j-1)*μ₁^(3n-2-i-j)*μ₂*μ₁^i*μ₂*μ₁^j*λ₁ for i=0:3n-2 for j=0:3n-2-i) + sum(μ₁^(3i)*μ₃*μ₁^(3n-1-3i)*λ₁ for i=0:n-1)
end

=#
