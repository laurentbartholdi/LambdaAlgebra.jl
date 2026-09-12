# julia --project=. benchmark/unstable_periodic.jl 100
# The argument bounds TOTAL degree, not topological degree.
using LambdaAlgebra
const LA=LambdaAlgebra
T=isempty(ARGS) ? 100 : parse(Int,ARGS[1])
T>=0 || error("total degree must be nonnegative")
A,_,_=periodic_algebra(p=3,dimension=3)
# Compile and warm up before measuring the larger diagonals.
LA.cache_basis!(A,min(T,20))
start=time_ns()
for n in unique(sort!([collect(40:20:T);T]))
    timing=@timed LA.cache_basis!(A,n)
    h=sum(count(LA.is_alive,b) for (d,b) in A.basis if d.top>0;init=0)
    println((total_degree=n,seconds=timing.time,cumulative_seconds=(time_ns()-start)/1e9,
             allocated_bytes=timing.bytes,positive_homology_classes=h,
             stats=LA.curtis_stats(A)))
    flush(stdout)
end
