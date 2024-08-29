using TeNe
include("../StatePreparation/gatesequence.jl")

## TFIM Parameters 
N = 24
d = 2
J = 0.4

# MPS/Circuit paramters 
D = 8
sub_layers = 2

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

### Decompose as staircase 
m = Int(ceil(log2(D))) + 1
ψ = productmps(Qubits(), ["dn" for _ = 1:N])
gates = randomstaircasecircuit(2, N, 1, m)
optim = StateOptimiser(ϕ, gates, ψ)

ψ′ = productmps(Qubits(), ["dn" for _ = 1:N])
applygates!(gates, ψ′)
println("Overlap Before = " * string(real(abs(inner(ψ′, ϕ))^2)))


# Decompose the staircase
U = Circuit(d, N) # The final circit 
for i = Base.OneTo(length(gates.layers))
    # Fetch the gate 
    g = gates.layers[i].gates[1]

    # Create subgates 
    if i == 1
        sub_gates = GateSequence(g, Tuple(1:m)...); # Each qubit meets a zero state 
    else
        sub_gates = GateSequence(g, m) # Only the final qubit meets a zero state
    end

    # Add the gauge
    if i != length(gates.layers)
        add!(sub_gates, TeNe._unitary_close_to_id(d, m-1), 2)
    end

    # Add 2-body gates 
    for k = 1:sub_layers
        for j = 1:m-1
            add!(sub_gates, TeNe._unitary_close_to_id(d, 2), j)
        end
    end

    # Do the optimisation 
    ols = optimise!(sub_gates)
    println("Fidelity of gate $(i): $(real(ols[end]))")

    # Absorb the gauge 
    if i != length(gates.layers)
        g′ = gates.layers[i+1].gates[1]
        new_gate = contract(
            tensor(g′),
            tensor(sub_gates.gatelist[begin]),
            Base.range(2, 2*m-2, step=2),
            Base.range(1, 2*m-3, step=2)
        )
        for j = 1:m-1
            new_gate = permutedim(new_gate, ndims(new_gate)-m+j+1, 2*j)
        end
        g′.gate .= new_gate
    end

    ### Add the gates to the full circuit 
    for j = Base.range(length(sub_gates.gatelist), i != length(gates.layers) ? 2 : 1, step=-1)
        add!(U, sub_gates.gatelist[j], Tuple(Base.range(i+sub_gates.qubits[j]-1, i+sub_gates.qubits[j])))
    end
end

ψ′ = productmps(Qubits(), ["dn" for _ = 1:N])
applygates!(U, ψ′)
println("Overlap After = " * string(real(abs(inner(ψ′, ϕ))^2)))