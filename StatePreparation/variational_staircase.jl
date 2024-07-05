using TeNe
using Printf

function variational_staircase(ϕ::GMPS, ψ::GMPS, m::Int; ϵ::Float64=1e-8, verbose::Bool=true)
    # Create a set of gates to optimise
    d = dim(ψ)
    gates = [TeNe._unitary_close_to_id(d, m, 0.1) for _ = 1:N-m+1]

    # Iterate until convergence 
    lastcost = _measure_overlap(ψ, ϕ, gates)
    costs = [lastcost]
    iter = 1
    while true
        for site = 1:N-m+1
            # Apply gates up to the one which is optimised.
            ψ′ = deepcopy(ψ)
            for i = 1:site-1
                applygate!(gates[i], ψ′, i; cutoff=1e-16)
            end
            
            # Apply gates after the one which is optimised.
            ϕ′ = deepcopy(ϕ)
            for i = 1:N-m+1-site
                gate = deepcopy(gates[end+1-i])
                sites = map(j->isodd(j) ? j + 1 : j - 1, 1:2*m)
                gate.gate = conj(permutedims(gate.gate, sites))
                applygate!(gate, ϕ′, N+1-i, true; cutoff=1e-16)
            end
    
    
            # Construct y 
            projϕ = ProjMPS(ϕ′, ψ′; center=site)
            prod = block(projϕ, site-1)
            rightblock = block(projϕ, site+m)
            for i = 1:m
                prod = contract(prod,  ϕ′[site+i-1], 2*(i-1)+1, 1, false, true)
                prod = contract(prod,  ψ′[site+i-1], 2*(i-1)+1, 1, false, false)
                prod = permutedim(prod, ndims(prod)-2, ndims(prod)-1)
            end
            prod = contract(prod, rightblock, (ndims(prod)-1, ndims(prod)), (1, 2))
            prod = conj(prod)
            new_gate = makeunitary(creategate(prod))
            gates[site] = new_gate
        end

        # Measure cost
        cost = _measure_overlap(ψ, ϕ, gates)
        push!(costs, cost)

        # Output 
        if verbose
            @printf("iter=%d, objective=%.8E, \n", iter, real(cost))
        end

        # Check for convergence
        if abs((cost - lastcost) / (0.5*(cost+lastcost))) < ϵ
            break
        end
        lastcost = cost
        iter += 1
    end
    
    return gates, costs
end

function _measure_overlap(ψ::MPS, ϕ::MPS, gates)
    ψ′ = deepcopy(ψ)
    for i = 1:length(gates)
        applygate!(gates[i], ψ′, i; cutoff=1e-16)
    end
    return abs(inner(ψ′, ϕ))^2
end

function _contract_gate(x, prod)
    prod = contract(prod, x, (1, 3), (2, 4), false, true)
    prod = permutedims(prod, (3, 1, 4, 2))
    return prod
end