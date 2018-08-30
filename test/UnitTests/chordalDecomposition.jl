using QOCS,SparseArrays,LinearAlgebra, FileIO

data = load("./ProblemFiles/decomp_testProblem.jld2")

P = data["P"]
q = data["q"]
A = data["A"]
b = data["b"]
objTrue = data["objTrue"]

cs = QOCS.Constraint(-A,b,QOCS.PositiveSemidefiniteCone())
model = QOCS.Model()
assemble!(model,P,q,[cs])
# # solve problem with sparsity turned off
settings = QOCS.Settings(decompose=false,obj_true=objTrue)
res = QOCS.optimize!(model,settings);


# # solve problem with sparsity exploitation
model = QOCS.Model()
assemble!(model,P,q,[cs])
settings_decomp = QOCS.Settings(decompose=true,obj_true=objTrue)
res_decomp = QOCS.optimize!(model,settings_decomp);


# # compare results
@testset "Chordal Decomposition" begin
  @test abs(res_decomp.objVal-res.objVal) < 1e-3
  @test abs(res_decomp.objVal-objTrue) < 1e-3
end