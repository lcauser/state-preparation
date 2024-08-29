using TeNe 

### System parameters; TFIM 
N = 30 # system size 
J = 0.4 # Iteraction


### Circuit optimisation 
χ = 8 # Determines how many qubits each gate in the staircase span 
m = Int(ceil(log2(χ))) + 1
sub_layers = 2 # Number of sub layers in the circuit 

### Do DMRG to find MPS 
# Create the Hamiltonain as an MPO
H = OpList(Qubits(), N)
for i = 1:N
    add!(H, "x", i, -1)
end
for i = 1:N-1
    add!(H, ["z", "z"], [i, i+1], -J)
end
H = MPO(H)

# Optimisation
ψ = randommps(2, N, 1)
energy, optim = dmrg(ψ, H; cutoff=1e-16, tol=1e-10)

### Create the circuit 
U = Circuit(2, N, CircuitMPS()) # The connector CircuitMPS makes it nice to work with MPS 
for i = 1:N+1-m
    for j = 1:sub_layers
        for k = 1:m-1
            add!(U, TeNe._unitary_close_to_id(2, 2, 1.0), (i+m-1-k, i+m-k)) # 1.0 is just the amount of noise
        end
    end
end

### Optimise the circuit to the MPS 
ϕ = productmps(Qubits(), ["dn" for _ = 1:N])
optim = StateOptimiser(ψ, U, ϕ)