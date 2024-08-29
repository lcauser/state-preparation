using TeNe 
import TeNe: add!

mutable struct GateSequence
    gate::TeNe.CircuitGate
    gatelist::Vector{TeNe.CircuitGate}
    qubits::Vector{Int}
    ends::Vector{Int}
end

"""
    GateSequence(gate::TeNe.CircuitGate)

Create a gate sequence that will be contracted with the given `gate`.
"""
function GateSequence(gate::TeNe.CircuitGate)
    return GateSequence(adjoint(gate), [], [], fill(0, length(gate)))
end
function GateSequence(gate::TeNe.CircuitGate, zeros::Int...)
    return GateSequence(adjoint(gate), [], [], Vector(map(j->(j in zeros) ? 2 : 0, Base.OneTo(length(gate)))))
end



"""
    add!(gates::GateSequence, gate::TeNe.CircuitGate, first::Int)

Add a gate to a sequence of gates. Use `first` to specify the first qubit which
the gate acts on.
"""
function add!(gates::GateSequence, gate::TeNe.CircuitGate, first::Int)
    if first + length(gate) - 1 <= length(gates.gate) && first > 0
        push!(gates.gatelist, gate)
        push!(gates.qubits, first)
    else
        throw(ArgumentError("The specified indices of the qubits exceed the length of the gate."))
    end
end


"""
    project(gates::GateSequence, which::Int)

Contract the overlap of the gate sequence with the specified gate, except for
at the index `which`.
"""
function project(gates::GateSequence, which::Int)
    # Fetch the gate 
    gate = gates.gate.gate

    # Contract the gates which appear before the specified gate
    for i = Base.OneTo(which-1)
        gate = _contract_after(gate, gates, i)
    end

    # Projections before 
    for i = Base.OneTo(length(gates.gate))
        if gates.ends[i] != 0
            gate = _project_before(gate, i, gates.ends[i])
        end
    end

    # Contract the gate which appear after the specified gate 
    for i = reverse(Base.range(which+1, length(gates.gatelist)))
        gate = _contract_before(gate, gates, i)
    end

    # Trace out the ends 
    for _ = Base.OneTo(gates.qubits[which]-1)
        gate = trace(gate, 1, 2)
    end
    offset = 2*length(gates.gatelist[which])
    for _ = Base.range(gates.qubits[which]+length(gates.gatelist[which]), length(gates.gate))
        gate = trace(gate, offset+1, offset+2)
    end

    return creategate(deepcopy(gate))
end

"""
    overlap(gates::GateSequence)

Calculate the overlap of sequence of gates with the specified gate.
"""
function overlap(gates::GateSequence)
    # Fetch the gate 
    gate = gates.gate.gate

    # Contract the gates which appear before the specified gate
    for i = Base.eachindex(gates.gatelist)
        gate = _contract_after(gate, gates, i)
    end

    # Trace out the ends
    for i = Base.OneTo(length(gates.gate))
        if gates.ends[i] == 0
            gate = trace(gate, 1, 2)
        else
            gate = gate[
                gates.ends[i], gates.ends[i],
                map(j->Base.range(Base.firstindex(gate, j), Base.lastindex(gate, j)),
                    Base.range(3, ndims(gate)))...
            ]
        end
    end
    return gate[] / 2^(sum(gates.ends .== 0))
end

function _contract_after(gate, gates, i)
    # Fetch the new gate 
    subgate = gates.gatelist[i]
    qubits = length(subgate)
    first_qubit = gates.qubits[i]

    # Contract with gate 
    gate = contract(
        gate, subgate.gate,
        Base.range(2*first_qubit, 2*(first_qubit-1+qubits), step=2),
        Base.range(1, 1+2*(qubits-1), step=2)
    )

    # Permute indices to the correct locations 
    for j = Base.OneTo(qubits)
        gate = permutedim(gate, ndims(gate)-qubits+j, 2*(first_qubit+j-1))
    end

    return gate
end

function _contract_before(gate, gates, i)
    # Fetch the new gate 
    subgate = gates.gatelist[i]
    qubits = length(subgate)
    first_qubit = gates.qubits[i]

    # Contract with the gate 
    gate = contract(
        subgate.gate,
        gate,
        Base.range(2, 2*qubits, step=2),
        Base.range(2*first_qubit-1, 2*(first_qubit-1+qubits)-1, step=2)
    )

    # Permute indices to the correct locations
    for j = reverse(Base.OneTo(qubits))
        gate = permutedim(gate, j, 2*(first_qubit-1+j)-1)
    end

    return gate
end

function _project_before(gate, i, val)
    ten = zeros(eltype(gate), size(gate, 2*i-1), size(gate, 2*i-1))
    ten[val, val] = 1
    gate = contract(gate, ten, 2*i-1, 1)
    gate = permutedim(gate, ndims(gate), 2*i-1)
    return gate
end

function optimise!(gates::GateSequence, ϵ::Float64=1e-8)
    # Measure convergence 
    lastcost = overlap(gates)
    overlaps = [lastcost]

    # Iterate
    while true
        # Update each gate
        for i in eachindex(gates.gatelist)
            new_gate = adjoint(makeunitary(project(gates, i)))
            gates.gatelist[i] = deepcopy(new_gate)
        end

        # Measure convergence
        cost = overlap(gates)
        push!(overlaps, cost)
        if abs((cost - lastcost) / (0.5*(cost+lastcost))) < ϵ
            break
        end
        lastcost = cost
    end
    
    return overlaps
end