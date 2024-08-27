using TeNe
include("../StatePreparation/variational_staircase.jl")
include("../StatePreparation/gatesequence.jl")

N = 30
D = 8
d = 2
J = 0.0

### Find the ground state to the Ising model using DMRG as an MPS 
H = OpList(Qubits(), N)
for i = 1:N
    add!(H, "x", i, -1)
end
for i = 1:N-1
    add!(H, ["z", "z"], [i, i+1], -J)
end
H = MPO(H)
ϕ = randommps(d, N, 1)
energy = dmrg(ϕ, H; cutoff=0, maxdim=D, tol=1e-10)

m = Int(ceil(log2(D))) + 1


qu = Qubits()
ψ = productmps(qu, ["up" for _ = 1:N])

gates, costs = variational_staircase(ϕ, ψ, m);

sub_gates_all = []
#for i = 1:length(gates)
    for i = 1:2
    # Create subgates 
    if i == 1
        sub_gates = GateSequence(gates[i], Tuple(1:length(gates[i]))...)
    else
        sub_gates = GateSequence(gates[i], length(gates[i]))
    end

    # Add the gauge
    if i != length(gates)
        add!(sub_gates, makeunitary(creategate(rand(ComplexF64, fill(2, 2*length(gates[i])-2)...))), 2)
    end

    # Add 2-body gates 
    for j = 1:length(gates[i])-1
        add!(sub_gates, makeunitary(creategate(rand(ComplexF64, 2, 2, 2, 2))), j)
    end

    if i != 1
        add!(sub_gates, makeunitary(creategate(rand(ComplexF64, fill(2, 2*length(gates[i])-2)...))), 1)
    end

    # Do the optimisation 
    ols = optimise!(sub_gates)
    println(real.(ols))

    # Absorb the gauge 
    if i != length(gates)
        new_gate = contract(
            tensor(gates[i+1]),
            tensor(sub_gates.gatelist[1]),
            Base.range(m+1, 2*m-1),
            Base.OneTo(m-1)
        )
        new_gate = permutedim(new_gate, m+1, 2*m)
        gates[i+1].gate = deepcopy(new_gate)
    end

    push!(sub_gates_all, sub_gates)
end