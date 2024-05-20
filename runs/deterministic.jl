using TeNe
using KrylovKit

function measure_overlap(ψ::MPS, ϕ::MPS, gates)
    ψ′ = deepcopy(ψ)
    for i = 1:N-1
        applygate!(gates[i], ψ′, i; cutoff=1e-16)
    end
    return abs(inner(ψ′, ϕ))^2
end

function contract_gate(x, prod)
    prod = contract(prod, x, (2, 4), (2, 4))
    prod = permutedims(prod, (1, 3, 2, 4))
    #prod = contract(prod, x, (2, 4), (2, 4))
    #prod = permutedims(prod, (3, 1, 4, 2))
    return prod
end

d = 2
N = 20 
D = 2

ϕ = randommps(d, N, D)
movecenter!(ϕ, 1)

gates = [TeNe._unitary_close_to_id(d, 2, 0.0) for _ = 1:N-1]

qu = Qubits()
ψ = productmps(qu, ["up" for _ = 1:N])

ψ′ = deepcopy(ψ)
for i = 1:N-1
    applygate!(gates[i], ψ′, i; cutoff=1e-16)
end
println("Overlap before: $(measure_overlap(ψ, ϕ, gates))")

ϕ′ = deepcopy(ϕ)
for i = 1:N-2
    applygate!(ϕ′, gates[end+1-i], N+1-i, true; cutoff=1e-16)
end


# Contruct M
projψ = ProjMPS(ψ, ψ)
leftblock = block(projψ, 0)
rightblock = block(projψ, 3)
prod = contract(leftblock, ψ[1], 1, 1, false, true)
prod = contract(prod, ψ[1], 1, 1, false, false)
prod = contract(prod, ψ[2], 2, 1, false, true)
prod = contract(prod, ψ[2], 3, 1, false, false)
M = contract(prod, rightblock, (4, 6), (1, 2); tocache=false)
f(x) = contract_gate(x, M)

# Construct y 
projϕ = ProjMPS(ψ, ϕ′)
leftblock = block(projϕ, 0)
rightblock = block(projϕ, 3)
prod = contract(leftblock, ψ[1], 1, 1, false, true)
prod = contract(prod, ϕ′[1], 1, 1, false, false)
prod = contract(prod, ψ[2], 2, 1, false, true)
prod = contract(prod, ϕ′[2], 3, 1, false, false)
y = contract(prod, rightblock, (4, 6), (1, 2); tocache=false)

new_gate, info = linsolve(f, y, gates[1].gate; ishermitian=false)