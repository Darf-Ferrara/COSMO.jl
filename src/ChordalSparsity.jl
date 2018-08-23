module ChordalSparsity


using ..QOCS,..Graphs, ..Trees, SparseArrays, LinearAlgebra
export SparsityPattern, findStackingMatrix, CliqueSet, getClique, chordalDecomposition!, reverseDecomposition!, findCommonSparsity
export vecToMatInd,findCommonSparsityPattern, findStackingMatrix


# ---------------------------
# STRUCTS
# ---------------------------
mutable struct SparsityPattern
  g::Graphs.Graph
  sntree::Trees.SuperNodeTree


  # constructor for sparsity pattern
  function SparsityPattern(A::Array{Int64,1},N::Int64,NONZERO_FLAG::Bool)
    g = Graph(A,N,NONZERO_FLAG)
    sntree = SuperNodeTree(g)

    return new(g,sntree)
  end

end



mutable struct CliqueSet
  cliqueInd::Vector{Int64} # an index set pointing to first element of each clique
  cliques::Vector{Int64} # just a vector with all cliques stacked
  vlen::Int64
  nBlk::Vector{Int64} # sizes of blocks
  N::Int64 # number of cliques in the set

  function CliqueSet(cliqueArr::Array{Array{Int64,1}})
    p = length(cliqueArr)
    # initialize fields
    cliqueInd = zeros(p)
    cliques = zeros(Int64,sum(map(x->length(x),cliqueArr)))
    nBlk = zeros(p)
    b = 1
    iii = 1
    for c in cliqueArr
      cliqueInd[iii] = b
      cliques[b:b+length(c)-1] = c
      nBlk[iii] = length(c)^2
      b +=length(c)
      iii +=1
    end

    return new(cliqueInd,cliques,length(cliques),nBlk,p)
  end
end



# ---------------------------
# FUNCTIONS
# ---------------------------
function _contains(convexSets::Array{AbstractConvexSet},t::Type{<:AbstractConvexSet})
  for set in convexSets
    if typeof(set) == t
      return true
    end
  end
  return false
end

function numSets(convexSets::Array{AbstractConvexSet},t::Type{<:AbstractConvexSet})
  number = 0
  for set in convexSets
    typeof(set) == t && (number+=1)
  end
  return number
end

function shiftIndices!(convexSets::Array{AbstractConvexSet},shift::Int64)
  for set in convexSets
      set.indices.start += shift
      set.indices.stop += shift
  end
end

function chordalDecomposition!(model::QOCS.Model,settings::QOCS.Settings,chordalInfo::QOCS.ChordalInfo)

  # do nothing if no psd cones present in the problem
  if _contains(model.convexSets,QOCS.PositiveSemidefiniteCone)
    settings.decompose = false
    return nothing
  end
  psdCones = filter(x->typeof(x) == QOCS.PositiveSemidefiniteCone,model.convexSets)

  # find sparsity pattern for each cone
  numCones = length(psdCones)
  spArr = Array{ChordalSparsity.SparsityPattern,1}(numCones)
  cliqueSets = Array{ChordalSparsity.CliqueSet,1}(numCones)

  # find sparsity pattern graphs and clique sets for each cone
  for (iii,cone) in enumerate(psdCones)
    ind = cone.indices
    csp = findCommonSparsity(A[ind,:],b[ind])
    cDim = Int(sqrt(cone.dim))
    sp = SparsityPattern(csp,cDim,true)
    spArr[iii] = sp
    cliqueSets[iii] = CliqueSet(sp.cliques)
  end

  # find transformation matrix H and store it
  H, decomposedPSDCones = findStackingMatrix(K,cliqueSets,model.m)
  chordalInfo.H = H

  # augment the system
  m,n = size(A)
  mH,nH = size(H)
  model.P = blkdiag(model.P,spzeros(nH,nH))
  model.q = vec([model.q;zeros(nH)])
  model.A = [model.A H; spzeros(nH,n) -sparse(1.0I,nH,nH)]
  model.b = vec([model.b;zeros(nH)])
  model.m = size(model.A,1)
  model.n = size(model.A,2)

  filter!(x->typeof(x) != QOCS.PositiveSemidefiniteCone,model.convexSets)
  model.convexSets = [model.convexSets;decomposedPSDCones]
  shiftIndices!(model.convexSets,mH)
  model.convexSets = [QOCS.Zeros(mH,1:mH);model.convexSets]
  nothing
end

# find the zero rows of a sparse matrix a
function zeroRows(a::SparseMatrixCSC,DROPZEROS_FLAG::Bool)
    DROPZEROS_FLAG && dropzeros!(a)
    passive = trues(a.m)
    for r in a.rowval
        passive[r] = false
    end
    return findall(passive)
end

function nzrows(a::SparseMatrixCSC,DROPZEROS_FLAG::Bool)
    DROPZEROS_FLAG && dropzeros!(a)
    active = falses(a.m)
    for r in a.rowval
        active[r] = true
    end
    return findall(active)
end

function numberOfOverlapsInRows(A::SparseMatrixCSC)
  # sum the entries row-wise
  numOverlaps = sum(A,2)
  ri = findall(x-> x > 1,numOverlaps)
  return ri, numOverlaps[ri]
end

function findCommonSparsity(A,b)
  AInd = ChordalSparsity.nzrows(A,false)
  # commonZeros = AInd[find(x->x==0,b[AInd])]
  bInd = findall(x->x!=0,b)
  commonNZeros = union(AInd,bInd)

  return commonNZeros
end

function findCommonSparsityPattern(Asub,bsub)
  m,n = size(Asub)
  AInd = zeroRows(Asub,false)
  commonZeros = AInd[find(x->x==0,b[AInd])]
  mSize = Int(sqrt(m))
  csp = spzeros(Int64,mSize,mSize)
  csp[:,:] = 1

  for ind in commonZeros
    i,j = vecToMatInd(ind,mSize)
    csp[i,j] = 0
  end
  return csp
end

function vecToMatInd(ind::Int64,n::Int64)
  ind > n^2 && error("Index ind out of range.")
  ind == 1 && (return 1,1)

  r = ind % n

  if r == 0
    j = Int(ind/n)
    i = n
  else
    j = Int(floor(ind/n) + 1)
    i = r
  end
  return i,j
end


function getClique(cs::ChordalSparsity.CliqueSet,ind::Int64)
  len = length(cs.cliqueInd)
  ind > len && error("Clique index ind=$(ind) is higher than number of cliques in the provided set:$(len).")
  ind < len ? (c = cs.cliques[cs.cliqueInd[ind]:cs.cliqueInd[ind+1]-1]) : (c = cs.cliques[cs.cliqueInd[ind]:end])
  return c
end

# function finds the transformation matrix H to decompose the vector s into its parts and stacks them into sbar, also returns  Ks
function findStackingMatrix(psdCones::Array{QOCS.PositiveSemidefiniteCone},cliqueSets::Array{ChordalSparsity.CliqueSet,1},m::Int64)

  numCones = length(psdCones)
  numCones != length(cliqueSets) && error("Number of psd cones and number of clique sets don't match.")

  stackedSizes = zeros(Int64,numCones)
  for iii=1:numCones
    stackedSizes[iii] = sum(cliqueSets[iii].nBlk)
  end

  bK = m - sum(map(x->x.dim,psdCones))

  # length of stacked vector sBar
  n = bK + sum(stackedSizes)

  H = spzeros(m,n)
  H[1:bK,1:bK] = speye(bK)
  bK += 1
  b = bK
  decomposedPSDCones = Array{QOCS.PositiveSemidefiniteCone}(undef,0)
  for kkk = 1:length(cliqueSets)
    cliqueSet = cliqueSets[kkk]
    mH = Int64
    for iii=1:cliqueSet.N
      mk = Int(sqrt(cliqueSet.nBlk[iii]))
      nk = Int(sqrt(psdCones[iii].dim))
      Ek = zeros(mk,nk)
      c = getClique(cliqueSet,iii)
      jjj = 1
      for v in c
        Ek[jjj,v] = 1
        jjj+=1
      end
      Hkt = (kron(Ek,Ek))'
      mH,nH = size(Hkt)
      H[b:b+mH-1,bK:bK+nH-1] = Hkt
      # create new cone and push to decomposedCone Array
      push!(decomposedPSDCones,QOCS.PositiveSemidefiniteCone(nH,bK:bK+nH-1))
      bK += nH
    end
    b += mH
  end
  return H, decomposedPSDCones

end

function reverseDecomposition!(ws::QOCS.Workspace,settings::QOCS.Settings)
  mO = ws.ci.originalM
  nO = ws.ci.originalN
  H = ws.ci.H

  sbar = ws.x[nO+1:end]
  s2 = ws.s[mO+1:end]

  ws.s = H*ws.s[mO+1:end]
  ws.x = ws.x[1:nO]

  # fill dual variables such that μ_k  = H_k μ for k=1,...,p
  fillDualVariables!(ws)
  # if user requests, perform positive semidefinite completion on entries of μ that were not in the decomposed blocks
  settings.completeDual && psdCompletion!(ws)
  δs = s2 - sbar
  maxRowH = maximum(sum(H,2))
  return δs, maxRowH
end

function fillDualVariables!(ws::QOCS.Workspace)
  mO = ws.ci.originalM
  H = ws.ci.H

  # this performs the operation μ = sum H_k^T *  μ_k causing an addition of (identical valued) overlapping blocks
  ws.μ = H*ws.μ[mO+1:end]
  ws.ν = H*ws.ν[mO+1:end]

  # to remove the overlaps we take the average of the values for each overlap by dividing by the number of blocks that overlap in a particular entry, i.e. number of 1s in each row of H
  rowInd,nnzs = numberOfOverlapsInRows(H)

  for iii=1:length(rowInd)
    ri = rowInd[iii]
    ws.μ[ri] = ws.μ[ri]/nnzs[iii]
    ws.ν[ri] = ws.ν[ri]/nnzs[iii]
  end
end

# complete the dual variable
function psdCompletion!(ws::QOCS.Workspace)
  return nothing
end


end #MODULE