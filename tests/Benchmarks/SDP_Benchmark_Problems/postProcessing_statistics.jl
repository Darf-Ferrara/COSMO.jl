workspace()
include("../../../src/Solver.jl")
include("../solverComparison/Compare.jl")
 include("./LatexExporter.jl")

# load data files
using Formatting, JLD, OSSDP, Compare
using LatexExport

folderName = "CompareSparseFormat_large_termination"
dir = "../../resultDataFiles/SDP_Benchmark_Problems/"*folderName
results = []
for f in filter(x -> endswith(x, ".jld"), readdir(dir))
    f = split(f,".")[1]
    push!(results,String(f))
end



  metricsArr = Array{Compare.ProblemMetrics}(length(results))

# Step 1: loop over all problem types and the combined file and calculate important metrics
for iii=1:length(results)
  r = results[iii]
  data = JLD.load(dir*"/"*r*".jld")
  resData = data["resData"]

  contains(resData[1].problemType, "Combine") ? COMBINE_FLAG = true : COMBINE_FLAG = false
  # create data object to hold the metrics for this problem type
  pm = Compare.ProblemMetrics(resData[1].problemType,COMBINE_FLAG,resData[1].ind,length(resData))

  # loop over solver and compute metrics
  k = 1
  for s in resData
    solvedInd = find(x->x==:solved,s.status)
    # calculate mean
    meanIterAll = mean(s.iter)
    length(solvedInd) > 0 ? meanIterSolved = mean(s.iter[solvedInd]) : meanIterSolved = Inf

    numSolved = length(solvedInd)
    percSolved = numSolved/s.ind

    meanErr = mean(abs.(s.objVal - s.objTrue)[solvedInd])
    meanRunTime = mean(s.runTime)
    meanIterTime = mean(s.iterTime)
    meanAvgIterTime = mean((s.iterTime./s.iter))*1000
    meanSetupTime = mean(s.setupTime)
    meanNZ = mean(s.problemDim[:,3])
    sm = Compare.SolverMetrics(s.solverName,s.adaptionON,s.scalingON,meanIterAll,meanIterSolved,numSolved,percSolved,meanErr,meanRunTime,meanIterTime,meanAvgIterTime,meanSetupTime,meanNZ)
    pm.solverResults[k] = sm
    k+=1
  end
  metricsArr[iii] = pm
end

# permute metricsArr to have combined results at the end of array
cind = find(x-> x.combinedData,metricsArr)
if length(cind) > 0
  cind = cind[1]
  if cind == 1
    p = [collect(2:length(metricsArr));1]
  elseif cind > 1 && cind < length(metricsArr)
    p = [collect(1:cind-1);collect(cind+1:length(metricsArr));2]
  end
  permute!(metricsArr,p)
end
# Step 2: Print results to screen
println("-"^80)
println("Postprocessing statistics:")
for pm in metricsArr
  println("-"^50)
  println(">>> Problem Type: $(pm.problemType)")
  println("- Combined Data: $(pm.combinedData)")
  println("- Number of Problems: $(pm.numProblems)")
  println("- Number of Solvers: $(pm.numSolvers)")

    println("Solver Name:\tNum Solved:\t% Solved:\tNonZ:\t\tIter (all):\tIter (solv):\tRT(all):\tIT(all):\tAvgIter[ms](all):\tSetupT(all):\tMean Error(solved):")
  for sm in pm.solverResults

    printfmt("{1:s}\t{2:d}\t\t{3:.2f}\t\t{4:.2f}\t\t{5:.2f}\t\t{6:.2f}\t\t{7:.4f}\t\t{8:.4f}\t\t{9:.4f}\t\t{10:.4f}\t\t{11:.4f}\n",sm.solverName[1:11],sm.numSolved,sm.percSolved,sm.meanNZ,sm.meanIterAll,sm.meanIterSolved,sm.meanRunTime,sm.meanIterTime,sm.meanAvgIterTime,sm.meanSetupTime,sm.meanErr)
    # println("Solver Name: $(sm.solverName[1:10])\tadaptionON: $(sm.adaptionON)\tscalingON: $(sm.scalingON)\tNum Solved: $(sm.numSolved)/$(pm.numProblems)\t% solved: $(sm.percSolved)\tMean Iter (all): $(sm.meanIterAll)\tMean Iter (solved): $(sm.meanIterSolved)")
    # println("Solver Name: $(sm.solverName)\t\tadaptionON: $(sm.adaptionON)\tscalingON: $(sm.scalingON)\tNum Solved: $(sm.numSolved)/$(pm.numProblems)\t% solved: $(sm.percSolved)\tMean Iter (all): $(sm.meanIterAll)\tMean Iter (solved): $(sm.meanIterSolved)")
  end
end
println("-"^80)

# Step 3: Print results to LaTeX table
resPath = dir*"/latex/"
!ispath(resPath) && mkdir(resPath)

 createLatexTable(metricsArr,resPath)
 # exportSolverSettings(metricsArr,resPath)



