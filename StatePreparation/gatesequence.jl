using TeNe 

mutable struct GateSequence
    gate::TeNe.CircuitGate
    gatelist::Vector{TeNe.CircuitGate}
    qubits::Vector{Int}
end

function GateSequence(gate::TeNe.CircuitGate)
    return GateSequence(adjoint(gate), [], [])
end

function add!(gates::GateSequence, gate::TeNe.CircuitGate, first::Int)
    if first + length(gate) - 1 <= length(gates.gate) && first > 0
        push!(gates.gatelist, gate)
        push!(gates.qubits, first)
    else
        throw(ArgumentError("The specified indices of the qubits exceed the length of the gate."))
    end
end


function project(gates::GateSequence, which::Int)
    # Create a copy of the gate


    # Do the gates which appear before the specified gate
    for i = Base.OneTo(which-1)

    end
end