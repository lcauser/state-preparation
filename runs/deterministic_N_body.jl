using TeNe
include("../StatePreparation/variational_staircase.jl")

d = 2
N = 30
D = 32

m = Int(ceil(log2(D))) + 1

ϕ = randommps(d, N, D)
movecenter!(ϕ, 1)

qu = Qubits()
ψ = productmps(qu, ["up" for _ = 1:N])

gates, costs = variational_staircase(ϕ, ψ, m);