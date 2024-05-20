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
    prod = contract(prod, x, (1, 3), (2, 4), false, true)
    prod = permutedims(prod, (3, 1, 4, 2))
    return prod
end

d = 2
N = 3
D = 2

ϕ = randommps(d, N, D)
movecenter!(ϕ, 1)

gates = [TeNe._unitary_close_to_id(d, 2, 0.1) for _ = 1:N-1]

qu = Qubits()
ψ = productmps(qu, ["up" for _ = 1:N])

for iter = 1:1
    for site = 2:2
        println("Overlap before: $(measure_overlap(ψ, ϕ, gates))")
        ψ′ = deepcopy(ψ)
        for i = 1:site-1
            applygate!(gates[i], ψ′, i; cutoff=1e-16)
        end

        ϕ′ = deepcopy(ϕ)
        for i = 1:N-1-site
            #applygate!(ϕ′, gates[end+1-i], N+1-i, true; cutoff=1e-16)
            gate = deepcopy(gates[end+1-i])
            gate.gate = conj(permutedims(gate.gate, (2, 1, 4, 3)))
            applygate!(gate, ϕ′, N+1-i, true; cutoff=1e-16)
        end


        # Contruct M
        projψ = ProjMPS(ψ′, ψ′; center=site)
        leftblock = block(projψ, site-1)
        rightblock = block(projψ, site+2)
        prod = contract(leftblock, ψ′[site], 1, 1, false, true)
        prod = contract(prod, ψ′[site], 1, 1, false, false)
        prod = contract(prod, ψ′[site+1], 2, 1, false, true)
        prod = contract(prod, ψ′[site+1], 3, 1, false, false)
        M = contract(prod, rightblock, (4, 6), (1, 2); tocache=false)
        f(x) = contract_gate(x, M)

        # Construct y 
        projϕ = ProjMPS(ϕ′, ψ′; center=site)
        leftblock = block(projϕ, site-1)
        rightblock = block(projϕ, site+2)
        prod = contract(leftblock, ϕ′[site], 1, 1, false, true)
        prod = contract(prod, ψ′[site], 1, 1, false, false)
        prod = contract(prod, ϕ′[site+1], 2, 1, false, true)
        prod = contract(prod, ψ′[site+1], 3, 1, false, false)
        y = contract(prod, rightblock, (4, 6), (1, 2); tocache=false)

        new_gate, info = linsolve(f, y, gates[site].gate; ishermitian=false, maxiter=1000)
        println(info)
        new_gate = makeunitary(creategate(new_gate))
        gates[site] = new_gate
        println("Overlap after: $(measure_overlap(ψ, ϕ, gates))")
    end
end