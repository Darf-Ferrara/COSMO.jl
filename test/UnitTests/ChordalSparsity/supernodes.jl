# Test Algorithm to check the supernode and clique algorithms
# Test the following properties:

# cliques
# - all cliques are maximal
# - the cliques cover all vertices
# - all cliques induce complete subgraphs


# snodes
# - first vertex v of each supernode is a representative vertex, i.e. every other vertex has a higher order
# - it holds snd(v) = {v, p(v), ..., p^n(v)} with deg+(v) = deg+(p^k(v)) + k, k=1,...,n, with n=|snd(v)|-1
# - the parent q(v) is the representative vertex of the supernode that contains the the parent of the highest element in snd(v)


using QOCS
using QOCS.Trees,QOCS.Graphs
using Test, Random, LinearAlgebra, SparseArrays

# create random seed (to get reproducable sequence of random numbers)
rng = MersenneTwister(123554);
# Define number of test matrices
nn = 100

function checkVertices(t,g::Graph)
  # all vertices of original graph have to occur in the supernodes of the super node elimination tree
  N = numberOfVertices(g)
  vertices = collect(1:N)
  numSN = length(t.snd_par)
  snptr = t.snptr
  for  iii=1:numSN
    if iii != numSN
      sn = t.snd[snptr[iii]:snptr[iii+1]-1]
    else
      sn = t.snd[snptr[iii]:end]
    end
    for v in sn
      deleteat!(vertices,findfirst(x->x == v,vertices))
    end
  end
  return size(vertices,1) == 0
end

function checkParents(t,g::Graph)
  # the parent of each supernode is the supernode that contains the parent of the highest order vertex
  N = length(t.snd_par)
  snd = t.snd
  snptr = t.snptr
  post = t.post
  postInv = Trees.invertOrder(post)
  par = t.par
  snd_par = t.snd_par
  for iii=1:N
    if iii != N
      sn = t.snd[snptr[iii]:snptr[iii+1]-1]
    else
      sn = t.snd[snptr[iii]:end]
    end

    highestVertex = sn[argmax(postInv[sn])]
    parentOfVertex = par[highestVertex]

    parInd = snd_par[iii]
    if parInd != 0
      if parInd != N
        snd_parent = t.snd[snptr[parInd]:snptr[parInd+1]-1]
      else
        snd_parent = t.snd[snptr[parInd]:end]
      end
      if !in(parentOfVertex,snd_parent)
        return false
      end
    end
  end
  return true
end

function repVertices(t)
  order = t.post
  orderInv = Trees.invertOrder(order)
  N = length(t.snptr)

  for iii=1:N
    snd = getSnd(t,iii)
    vRep = snd[1]

    for o_v in snd
      if o_v != vRep
        if orderInv[o_v] <= orderInv[vRep]
          return false
        end
      end
    end
  end
  return true
end

function sndDegrees(t,g)
  degrees = Graphs.higherDegrees(g)
  N = length(t.snptr)
  par = t.par

  for iii=1:N
    snd = getSnd(t,iii)
    vRep = snd[1]
    par_v = vRep
    for jjj=1:length(snd)
      if jjj>=2
        par_v = par[par_v]
        # #check deg+(v) = deg+(p^k(v)) + k, k=1,...,n, with n=|snd(v)|-1
        if !(degrees[vRep] == degrees[par_v] + (jjj-1))
          return false
        end
      end
    end
  end
  return true
end

# all vertices in a clique have to be pair-wise adjacent, i.e. induce a complete subgraph
function completeSubgraphs(t,g::Graph)
  Nc = length(t.chptr)
  chptr = t.chptr
  for iii = 1:Nc
    if iii != Nc
      clique = t.cliques[chptr[iii]:chptr[iii+1]-1]
    else
      clique = t.cliques[chptr[iii]:end]
    end

    if size(clique,1) > 1
      for v in clique
        otherNodes = filter(x->x!=v,clique)
        for cliqueNode in otherNodes
          if !in(cliqueNode,g.adjacencyList[v])
            println("iii:$(iii)")
            return false
          end
        end
      end
    end
  end
  return true
end

# all cliques are maximal, i.e. no clique contains another clique
function maximalCliques(t)
  cliques = t.cliques
  chptr = t.chptr
  Nc = length(chptr)
  for iii=1:Nc
    ch = getClique(t,iii)

    for jjj = 1:Nc
      if jjj != iii
        o_ch = getClique(t,jjj)
        if issubset(ch,o_ch)
          return false
        end
      end
    end
  end

  return true

end

@testset "SuperNode Tree Derivation from random Matrices" begin

for lll=1:nn

  # take random dimension
  #dim = rand(rng,20:100)
  dim = 8
  density = rand(rng,0.1:0.1:0.6)
  # create random sparse matrix
  A = sprand(rng,dim,dim,density)
  A = A+A'

  # create graph from A and make Graph chordal
  g = Graph(A)

  # create Tree from Graph
  sntree = Trees.SuperNodeTree(g)


  @test checkVertices(sntree,g)
  @test checkParents(sntree,g)
   @test repVertices(sntree)
   @test sndDegrees(sntree,g)
   @test completeSubgraphs(sntree,g)
  @test maximalCliques(sntree)
end
end




