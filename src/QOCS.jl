#__precompile__()
module QOCS

#export MathOptInterfaceQOCS
export optimize!, reset!, assemble!, warmStart!, Settings, Model, Result, Constraint, AbstractConvexSet
using SparseArrays,LinearAlgebra

include("./Algebra.jl")
include("./Helper.jl")
include("./Sets.jl")
include("./Types.jl")
include("./KKT.jl")
include("./Scaling.jl")
include("./Projections.jl")
include("./Residuals.jl")
include("./Parameters.jl")
include("./Infeasibility.jl")
include("./Printing.jl")
include("./Setup.jl")

using .Projections, .Scaling, .Parameters, .Infeasibility, .Residuals, .Printing, .Setup
include("./Solver.jl")
include("./Interface.jl")

end
