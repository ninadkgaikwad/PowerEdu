# Optimal Power Flow Functions

using Symbolics

include("Helper_Functions.jl")


function solveForEconomicDispatch(CDF_DF_File_pu::Vector{DataFrame},
    x::Vector{Num}, f::Num, h::Vector{Num};
    verbose::Bool=false)

    grad_f = Symbolics.gradient(f, x)
    jacobian_h = Symbolics.jacobian(h, x)
    @variables lambda[1:length(h)]  # Lagrange multipliers
    if length(h) == 1
        myprintln(verbose, "Okay, just a single equality.")
        lambda = lambda[1];
        jacobian_h = vec(jacobian_h)
    end
    # @variables lambda
    # @show jacobian_h
    KKT_conditions = vec([grad_f - jacobian_h*lambda; h])
    myprintln(verbose, KKT_conditions)
    # println([x; lambda])
    # println(KKT_conditions)
    solutions = Symbolics.solve_for(KKT_conditions, [x; lambda])

    myprintln(true, "Generator 1: $(solutions[1]) MW")
    myprintln(true, "Generator 2: $(solutions[2]) MW")
    myprintln(true, "Marginal Cost of Generation = \$ $(solutions[3]).")
    return solutions
end

# PL = 259;
# @variables P1 P2;
# x = [P1, P2];
# f₁ = 0.004P1^2 + 8P1;
# f₂ = 0.0048P2^2 + 6.4P2;
# h₁ = P1 + P2 - Pₗ₁;
# h₂ = h₁;
# # h = [h₁];
# h = [h₁, h₂];
# f = f₁ + f₂;
# solutions1 = solveForEconomicDispatch(dfpu, x, f, h, verbose=true);
# P1′, P2′, lambda₁′ = solutions1;
