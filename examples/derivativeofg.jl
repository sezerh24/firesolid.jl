
using Symbolics

# Define the symbolic variable
@variables ξ

# Define the constant ζ
ζ = 20

W = 1.78
# Define the function
g(ξ) = W * ξ^2 * (1 - ξ)^2/16

# Compute the derivative with respect to ξ
dg_dξ = Symbolics.derivative(g(ξ), ξ)

# Convert the derivative to a callable function
dg_dξ_func = Symbolics.build_function(dg_dξ, ξ)


# Define the callable function
dg_dξ_callable = eval(dg_dξ_func)

ξ = sort(rand(1000));

using Plots
plot(ξ,g.(ξ))
plot!(ξ,dg_dξ_callable.(ξ))


