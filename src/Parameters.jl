module Parameters
using ..QOCS,..Residuals, ..KKT
export setRhoVec!, adaptRhoVec!, updateRhoVec!


# set initial values of rhoVec
  function setRhoVec!(ws::QOCS.Workspace,settings::QOCS.Settings)
    p = ws.p
    #nEQ = sum(map(x->x.dim,filter(x->typeof(x) == QOCS.Zeros(), ws.p.convexSets))) #p.K.f
    #nINEQ = p.m - nEQ
    ws.ρ = settings.rho
    ws.ρVec = ws.ρ*ones(p.m)
    for set in ws.p.convexSets
      if typeof(set) == QOCS.Zeros
        ws.ρVec[set.indices] *= 1000
      end
    end
    #ws.ρVec = [1e3*ws.ρ*ones(nEQ);ws.ρ*ones(nINEQ)] #ws.ρ*ones(p.m)
    push!(ws.Info.rho_updates,ws.ρ)
    return nothing
  end


  # adapt rhoVec based on residual ratio
  function adaptRhoVec!(ws::QOCS.Workspace,settings::QOCS.Settings)
    # compute normalized residuals based on the working variables (dont unscale)
    IGNORE_SCALING = true
    r_prim, r_dual = calculateResiduals(ws,settings,IGNORE_SCALING)
    maxNormPrim, maxNormDual = maxResComponentNorm(ws,settings,IGNORE_SCALING)
    r_prim = r_prim/(maxNormPrim + 1e-10)
    r_dual = r_dual/(maxNormDual + 1e-10)

    newRho = ws.ρ * sqrt(r_prim/(r_dual+1e-10))
    newRho = minimum([maximum([newRho,settings.RHO_MIN]),settings.RHO_MAX])
    # only update rho if significantly different than current rho
    # FIXME: Should it be settings.rho or previous rho???
    if (newRho > settings.adaptive_rho_tolerance*ws.ρ) || (newRho < (1 ./ settings.adaptive_rho_tolerance)*ws.ρ)
      updateRhoVec!(newRho,ws,settings)
    end
    return nothing
  end

  function updateRhoVec!(newRho::Float64,ws::QOCS.Workspace,settings::QOCS.Settings)
    p = ws.p
    #nEQ = sum(map(x->x.dim,filter(x->typeof(x) == QOCS.Zeros(), ws.p.convexSets))) #p.K.f
    #nINEQ = p.m - nEQ

    ws.ρ = newRho
    ws.ρVec = ws.ρ*ones(p.m)
        for set in ws.p.convexSets
          if typeof(set) == QOCS.Zeros()
            ws.ρVec[set.indices] *= 1000
          end
        end    # log rho updates to info variable
    push!(ws.Info.rho_updates,newRho)
    factorKKT!(ws,settings)
    return nothing
  end

end #MODULE