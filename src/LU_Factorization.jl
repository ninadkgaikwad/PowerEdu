# LU_Factorization.jl

"""
    Compute_PLU(A, Tol_Num)

Computes LU decomposition using partial pivoting.

'''
# Arguments
- 'A': A square matrix
- 'Tol_Num': Tolerance for being near zero
'''
'''
# Output
- 'P': Permutation Matrix which keeps tabs on row ordering.
- 'L': Lower Triangular matrix of the LU Decomposition.
- 'U': Upper Triangular matrix of the LU Decomposition.
'''
"""
function Compute_PLU(A, Tol_Num)

    # Get size of A
    RowNum_A, ColNum_A = size(A)

    # Initializing P, L and U
    U = A
    L = Matrix{Float64}(I, RowNum_A, ColNum_A)
    P = Matrix{Float64}(I, RowNum_A, ColNum_A)

    # Computing P, L, U
    for ii in 1:RowNum_A # Through all the rows

        # Permuting Rows for a Zero valued Pivot
        for jj in ii+1:RowNum_A # Through all the rows

            # Checking if Pivot is near Zero
            if (!isapprox(U[ii,ii],0.0; atol=Tol_Num))

                break

            end

            # Storing ii Row in Temp Variables
            Temp_Row_U_ii = U[ii,:]
            Temp_Row_P_ii = P[ii,:]

            # Updating ii Row with jj Row
            U[ii,:] = U[jj,:]
            P[ii,:] = P[jj,:]

            # Updating jj Row with ii Row from Temp Variable
            U[jj,:] = Temp_Row_U_ii
            P[jj,:] = Temp_Row_P_ii

        end

        # Computing Current Pivot Factors
        CurrentPivot_Factor = U[ii+1:end,ii]/U[ii,ii]

        # Updating L using Current Pivot Factors
        L[ii+1:end,ii] = CurrentPivot_Factor

        # Updating U using Current Pivot Factors
        U[ii+1:end,:] -= CurrentPivot_Factor.*U[ii,:]'

    end

    return P, L, U

end

"""
    LU_ForwardSubstitution(L, b)

Computes y from L and b in the LU decomposition scheme.

'''
# Arguments
- 'L': Lower Triangular Matrix of LU Decomposition
- 'b': Vector of Knowns from the original problem
'''
'''
# Output
- 'y': Solution vector to the problem Ly = b.
'''
"""
function LU_ForwardSubstitution(L, b)

    # Getting Shape of L
    Row_L, Col_L = size(L)

    # Initializing: Vector y
    y = zeros(Row_L,1)

    # Initializing: Forward Substitution
    y[0] = b[0]/L[0,0]

    # Solving: Remaining Forward Substitution
    for ii in 2:Row_L

        y[ii] = (b[ii] - L[ii,1:ii-1]'*y[1:ii-1])/L[ii,ii]

    end

    return y

end

"""
    LU_BackwardSubstitution(U, y)

Computes x from U and y in the LU decomposition scheme.

'''
# Arguments
- 'U': Upper Triangular Matrix of LU Decomposition
- 'y': Solution vector to the problem Ly = b obtained from forward substitution
'''
'''
# Output
- 'x': Solution vector to the original problem Ax = b.
'''
"""
function LU_BackwardSubstitution(U, y)

    # Getting Shape of U
    Row_U, Col_U = size(U)

    # Initializing: Vector x
    x = zeros(Row_U,1)

    # Initializing: Backward Substitution
    x[end] = y[end]/U[end,end]

    # Solving: Remaining Backward Substitution
    for ii in Row_U-1:-1:1

        x[ii] = (y[ii] - U[ii,ii+1:end]'*x[ii+1:end])/U[ii,ii]

    end

    return x

end

"""
    PLU_Solve(A, Tol_Num)

Computes solution of Ax = b using LU decomposition by partial pivoting.

'''
# Arguments
- 'A': A square matrix
- 'b': Known vector
- 'Tol_Num': Tolerance for being near zero
'''
'''
# Output
- 'x': Solution vector to the original problem Ax = b.
'''
"""
function PLU_Solve(A, b, Tol_Num)

    # Computing the Permutation and LU Decomposition
    P, L, U = Compute_PLU(A, Tol_Num)

    # Effect of P on b
    b_Permuted = P*b

    # Computing y using Forward Susbtitution
    y = LU_ForwardSubstitution(L, b_Permuted)

    # Computing x using Backward Sustition
    x_Permuted = LU_BackwardSubstitution(U, y)

    # Unpermuting x_Permuted
    x = zeros(length(x_Permuted),1)

    for ii in 1:length(x_Permuted)

        # Finding Position of 1 in ii Row of P
        OriginalIndex = findall(x -> x == 1, P[ii,:])

        # Updating x from x_Permuted using OriginalIndex
        x[OriginalIndex[1],1] = x_Permuted[ii,1]

    end

    return x

end
