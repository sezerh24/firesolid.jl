using Symbolics

# Define the symbolic variable
@variables ξ

# Define the constant ζ
ζ = 20

# Define the function
h(ξ) = exp(ζ * (ξ - 1/2)) / (1 + exp(ζ * (ξ - 1/2)))

# Compute the derivative with respect to ξ
dh_dξ = Symbolics.derivative(h(ξ), ξ)

# Convert the derivative to a callable function
dh_dξ_func = Symbolics.build_function(dh_dξ, ξ)

# Print the generated function code
println(dh_dξ_func)

# Define the callable function
dh_dξ_callable = eval(dh_dξ_func)

# Test the function
ξ_val = 15.0
println("The derivative at ξ = $ξ_val is ", dh_dξ_callable(ξ_val))

ξ_val = sort(rand(1000));

using Plots
plot(ξ_val,h.(ξ_val))
plot!(ξ_val,dh_dξ_callable.(ξ_val))


