include("../StatePreparation/gatesequence.jl")

U = makeunitary(creategate(rand(ComplexF64, 2, 2, 2, 2, 2, 2, 2, 2)))

#U = makeunitary(creategate(rand(ComplexF64, 2, 2, 2, 2, 2, 2, 2, 2)))
#gates_test = GateSequence(U)
#add!(gates_test, U, 1)
#println(overlap(gates_test))

gates = GateSequence(U)
add!(gates, makeunitary(creategate(rand(ComplexF64, 2, 2, 2, 2, 2, 2))), 2)
add!(gates, makeunitary(creategate(rand(ComplexF64, 2, 2, 2, 2))), 1)
add!(gates, makeunitary(creategate(rand(ComplexF64, 2, 2, 2, 2))), 2)
add!(gates, makeunitary(creategate(rand(ComplexF64, 2, 2, 2, 2))), 3)
add!(gates, makeunitary(creategate(rand(ComplexF64, 2, 2, 2, 2))), 1)
add!(gates, makeunitary(creategate(rand(ComplexF64, 2, 2, 2, 2))), 2)
add!(gates, makeunitary(creategate(rand(ComplexF64, 2, 2, 2, 2))), 3)

costs = optimise!(gates)