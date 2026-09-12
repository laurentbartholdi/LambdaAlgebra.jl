# julia --project=. benchmark/compare_variants.jl 40 60 80 100
# End-to-end construction and homology, after warming both implementations.
# Cutoffs bound total degree mu+lambda+top; all mu degrees are included.
using LambdaAlgebra
const LA=LambdaAlgebra

function compute(constructor,T)
    A=constructor(p=3,dimension=3,top_degree=T)[1]
    A,homology(A)
end

function measure_periodic(T)
    GC.gc()
    r=@timed compute(periodic_algebra,T)
    A,H=r.value
    signature=Dict(d=>length(b) for (d,b) in H.basis if d.top>0 && !isempty(b))
    stats=LA.curtis_stats(A) # periodic only; ordinary uses no context engine
    (seconds=r.time,allocated_bytes=r.bytes,gc_seconds=r.gctime,
     public_entries=sum(length,values(A.basis)),contexts=stats.contexts,
     context_entries=stats.stored,words=stats.words,word_slots=stats.word_slots,
     positive_classes=sum(values(signature))),signature
end

function measure_usual(T)
    GC.gc()
    r=@timed compute(lambda_algebra,T)
    A,H=r.value
    signature=Dict(d=>length(b) for (d,b) in H.basis if d.top>0 && !isempty(b))
    (seconds=r.time,allocated_bytes=r.bytes,gc_seconds=r.gctime,
     public_entries=sum(length,values(A.basis)),positive_classes=sum(values(signature))),signature
end

compute(lambda_algebra,24)
compute(periodic_algebra,24)
cutoffs=isempty(ARGS) ? [40,60,80,100] : parse.(Int,ARGS)
for T in cutoffs
    T>=0 || error("total degree must be nonnegative")
    ordinary,hs=measure_usual(T)
    println((total_degree=T,variant="lambda-mu",ordinary...))
    flush(stdout)
    periodic,hp=measure_periodic(T)
    println((total_degree=T,variant="periodic",periodic...))
    @assert hp==hs "Positive-degree homology signatures differ at cutoff $T"
    println((total_degree=T,matching_positive_homology=true))
    flush(stdout)
end
