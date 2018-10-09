module Printing
  using ..QOCS,Printf, LinearAlgebra
  export printHeader, printResult, printIteration

  function printHeader(ws::QOCS.Workspace,settings::QOCS.Settings)
    n = ws.p.n
    m = ws.p.m
    K = ws.p.K
    numSOC = numSOCBlk = numSDP = numSDPBlk = 0
    if length(K.q) > 0
      numSOC = sum(K.q)
      numSOCBlk = length(K.q)
    else
      numSOC = 0
      numSOCBlk = 0
    end
    if length(K.s) > 0
      numSDP = sum(K.s)
      numSDPBlk = length(K.s)
    else
      numSDP = 0
      numSDPBlk = 0
    end

    settings.scaling > 0 ? scalingStatus = "on" : scalingStatus = "off"
    nnzInP = count(!iszero,ws.p.P) - count(!iszero,diag(ws.p.P)) + n
    nnzInM = 2*count(!iszero,ws.p.A) + nnzInP + m
    println("-"^66 * "\n" * " "^13 * "Quadratic Objective Conic Solver (QOCS)\n" * " "^25 * "Michael Garstka\n"  * " "^16 * "University of Oxford, 2017 - 2018\n" * "-"^66 * "\n")
    println("Problem:  x ∈ R^{$(n)},\n          constraints: A ∈ R^{$(m)x$(n)} ($(count(!iszero,ws.p.A)) nnz), b ∈ R^{$(m)},\n          matrix size to factor: $(n+m)x$(n+m) ($((n+m)^2) elem, $(nnzInM) nnz)")
    for (iii,set) in enumerate(sort(ws.p.convexSets,by=x -> -x.dim))
      setName = split(string(typeof(set)),".")[end]
      iii == 1 ? println("Sets:"*" "^5*"$(setName) of dim: $(set.dim)") : println(" "^10*"$(setName) of dim: $(set.dim)")
      if iii > 5
        print(" "^10*"... and $(length(ws.p.convexSets)-5) more")
        break
      end
    end
    println("Settings: ϵ_abs = $(@sprintf("%.1e",settings.eps_abs)), ϵ_rel = $(@sprintf("%.1e",settings.eps_rel)),\n" * " "^10 * "ϵ_prim_inf = $(@sprintf("%.1e",settings.eps_prim_inf)), ϵ_dual_inf = $(@sprintf("%.1e",settings.eps_dual_inf)),\n" * " "^10 * "ρ = $(settings.rho), σ = $(settings.sigma), α = $(settings.alpha),\n" * " "^10 * "max_iter = $(settings.max_iter),\n" * " "^10 * "scaling iter = $(settings.scaling) ($(scalingStatus)),\n" * " "^10 * "check termination every $(settings.check_termination) iter,\n" * " "^10 * "check infeasibility every $(settings.check_infeasibility) iter")
    println("Setup Time: $(round.(ws.times.setupTime*1000;digits=2))ms\n")
    println("Iter:\tObjective:\tPrimal Res:\tDual Res:\tRho:")
    nothing
  end

  function printIteration(ws::QOCS.Workspace,settings::QOCS.Settings,iter::Int64,cost::Float64,r_prim::Float64,r_dual::Float64)

    if mod(iter,1) == 0 || iter == 1 || iter == 2 || iter == settings.max_iter
      if mod(iter,settings.check_termination) == 0
        @printf("%d\t%.4e\t%.4e\t%.4e\t%.4e\n", iter,cost,r_prim,r_dual,ws.ρ)
      else
        @printf("%d\t%.4e\t ---\t\t\t---\n", iter,cost)
      end
    end
    nothing
  end


  function printResult(status::Symbol,iter::Int64,cost::Float64,rt::Float64)
    println("\n" * "-"^66 * "\n>>> Results\nStatus: $(status)\nIterations: $(iter)\nOptimal objective: $(round.(cost;digits=4))\nRuntime: $(round.(rt;digits=3))s ($(round.(rt*1000;digits=2))ms)\n")
    nothing
  end
end