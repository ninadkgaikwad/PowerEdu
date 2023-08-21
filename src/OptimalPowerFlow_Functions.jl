# Optimal Power Flow Functions

using Symbolics
using ForwardDiff
using NLsolve

function objective(p)
    P1, P2 = p
    F1 = 0.004 * P1^2 + 8 * P1
    F2 = 0.0048 * P2^2 + 6.4 * P2
    return F1 + F2
end

gradient_objective = p -> ForwardDiff.gradient(objective, p)

function system_of_equations!(F, p)
    grad = gradient_objective(p)
    F[1] = grad[1] - grad[2]
    F[2] = p[1] + p[2] - 259.0
end

initial_guess = [130.0, 130.0]  # starting with a simple guess
result = nlsolve(system_of_equations!, initial_guess);
optimal_p = result.zero

# function solveForEconomicDispatch(CDF_DF_File_pu::Vector{DataFrame},
#     x:Vector, f;
#     verbose=false)

#     df_dx = Symbolics.jacobian([f], x);
#     solutions = Symbolics.solve_for(df_dx, x);
#     return solutions
# end

function objective(p)
    f[1] = 0.004p[1] + 8p[1]
    f[2] = 0.0048p[2] + 6.4p[2]
    f = sum(f)
    return f
end

