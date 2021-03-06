# Helper Script that reads in the problems in the SDPA sparese file format and recreates and
# stores the problem data as matrices in Julia compatible JLD format

# The following is stored in the JLD file:
# m: dimension of x
# n: dimension of Fs,X
# nblocks: number of blocks in the diagonal structure of the matrices
# blockvec: vector of numbers that give the sizes on the individual blocks
# c: objective vector
# F: array containing the constraint matrices, F[1]=F0, F[2]=F1, F[m+1]=Fm
# optval: Optimal objective value (extracted from README file)

# by Michael Garstka
# University of Oxford, Control Group
# January 2017


workspace()
using JLD
# Specify path to .dat-s files
dirPath = "./sdplib/"

# create array of all relevant filenames that contain problems
fileNames = []
for f in filter(x -> endswith(x, "dat-s"), readdir(dirPath))
    f = split(f,".")[1]
    push!(fileNames,f)
end

#fileNames = ["qap10"]
for file in fileNames
  # open file and read everything into memory
  f = open(dirPath*file*".dat-s");
  lines = readlines(f)
  close(f)


  # initialize variables
  counter = 1
  # size of x
  m = Inf
  # size of symmetric F matrices
  n = Inf
  nblocks = Inf
  blockvec = zeros(1)
  c = zeros(1)
  F = zeros(1)
  currentM = -1
  Fm = zeros(1)


  for ln in lines
    # dont do anything if line contains a comment
    if ln[1] == '"' || ln[1] == "*"
      counter == 1
    else
      # m: number of constraint matrices (in SDPA it starts from 0 -> +1)
      if counter == 1
        m = parse(Int64,strip(ln)) + 1
        # now that m is known, initalize the array that holds all the constraint matrices F
        F = Array{SparseMatrixCSC}(m)

      # nblocks: number of blocks in the diagonal structure of the matrices
      elseif counter == 2
        nblocks = parse(Int64,strip(ln))

      # vector of numbers that give the sizes on the individual blocks
      # negative number indicates diagonal submatrix
      elseif counter == 3
      bvec = split(replace(strip(ln,['{','}','(',')']),"+",""))
       blockvec = [parse(Float64, ss) for ss in bvec]
       n = Int(sum(abs.(blockvec)))

       # objective function vector c
      elseif counter == 4
        # try first to split by whitespace as delimiter (also remove certain trailing, ending characters)
        cvec = split(replace(strip(ln,['{','}','(',')']),"+",""))
        # otherwise try comma as delimiter
        if length(cvec) == 1
          cvec = split(replace(strip(ln,['{','}','(',')']),"+",""),",")
        end
        c = [parse(Float64, ss) for ss in cvec]

      # all other lines contain constraint matrices with one entry per line
      # save them directly as sparse matrix
      else
        # FIXME: Accuracy
        line = [parse(Float64, ss) for ss in split(ln)]
        matno = Int(line[1])
        blkno = Int(line[2])
        i = Int(line[3])
        j = Int(line[4])
        entry = line[5]

        if matno > currentM
          # save previously generated matrix in sparse format
          if currentM >= 0
            F[currentM+1] = sparse(Fm)
          end
          # create new temporary matrix that is filled
          Fm = zeros(n,n)
          currentM = matno
        end

        if blkno == 1
          iF = i
          jF = j
        else
          shift = Int(sum(abs.(blockvec[1:blkno-1])))
          iF = i+shift
          jF = j+shift
        end
        Fm[iF,jF] = entry
        # ensure symmetry
        if iF != jF
          Fm[jF,iF] = entry
        end
      end
      counter += 1
    end
  end
  # after last line has reached save the last F matrix
  F[currentM+1] = sparse(Fm)

  # extract solution for objective value from README file
  f = open(dirPath*"README")
  lines = readlines(f)
  close(f)
  lines = lines[20:113]
  ln = filter(s -> contains(s,file) ,lines)
  if length(ln) == 0
    optVal = "Not provided"
  else
    str = split(ln[1])[4]
    if contains(str,"primal")
      optVal = Inf
    elseif contains(str, "dual")
      optVal = -Inf
    else
    optVal = parse(Float64, str)
    end
  end

  # save to JLD file
  JLD.save(dirPath*file*".jld", "m", m-1,"n",n,"nblocks",nblocks,"blockvec",blockvec,"c",c,"F",F,"optVal",optVal)
  println("Saved problem: $(file) to "*dirPath*file*".jld, Optimal Value: $(optVal)")
end