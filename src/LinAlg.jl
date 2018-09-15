module LinAlg

  using LinearAlgebra, SparseArrays

  export left_right_diag_mat!, left_diag_mat!, right_diag_mat!,left_right_diag_sparse!

  # implements A = D1 * A * D2 where d1, d2 are vectors of the diagonal entries of D1, D2
  function left_right_diag_mat!(A::Array{Float64,2},d1,d2)
    # loop through cols
    for j = 1:size(A,2)
      # loop through rows
      for i = 1:size(A,1)
        A[i,j] = d1[i] * A[i,j] * d2[j]
      end
    end
  end

  function left_right_diag_sparse!(A::SparseMatrixCSC{Float64,Int64},d1,d2)
  # loop through cols
  for j = 1:size(A,2)
    b = A.colptr[j]
    e = A.colptr[j+1] - 1
    for i = b:e
      A.nzval[i] =  d1[A.rowval[i]] * A.nzval[i] * d2[j]
    end
  end
end

    # implements A = D * A where d, is a vector of the diagonal entries of D
  function left_diag_mat!(A::SparseMatrixCSC{Float64,Int64},d)
    # loop through cols
    for j = 1:size(A,2)
      # loop through rows
      for i = 1:size(A,1)
        A[i,j] = d[i] * A[i,j]
      end
    end
  end

      # implements A = A * D where d, is a vector of the diagonal entries of D
  function right_diag_mat!(A::SparseMatrixCSC{Float64,Int64},d)
    # loop through cols
    for j = 1:size(A,2)
      # loop through rows
      for i = 1:size(A,1)
        A[i,j] =  A[i,j] * d[j]
      end
    end
  end

end