
# Declare the type of the NamedTuple
NamedTupleType = NamedTuple{(:zama, :amam), Tuple{Function, Float64}}

# Define a function
func = (α, β) -> α * β

# Create a NamedTuple with the declared type
data = NamedTupleType((zama = func, amam = 42.0))

# Access the elements of the NamedTuple
println(data.zama)  # Prints the function
println(data.amam)  # Prints 42.0


# Define a function to handle a NamedTuple with specific fields
function process(data::NamedTuple{(:zama, :amam), Tuple{Function, Float64}})
    println("Processing NamedTuple with function and Float64")
    println("Function output: ", data.zama(2, 3))  # Example function call
    println("Float64 value: ", data.amam)
end

# Define a function to handle a NamedTuple with different fields
function process(data::NamedTuple{(:zama, :amam, :extra), Tuple{Function, Float64, Int}})
    println("Processing NamedTuple with function, Float64, and Int")
    println("Function output: ", data.zama(2, 3))  # Example function call
    println("Float64 value: ", data.amam)
    println("Extra Int value: ", data.extra)
end

# Define a function to handle other NamedTuples
function process(data::NamedTuple)
    println("Processing a general NamedTuple")
    println("Contents: ", data)
end


 # Species in the gas phase 
gasspec   = ["α", "β", "γ", "θ", "ζ", "η", "ϕ"]

#species in solid phase
solidspec = ["ξ", "μ"]
spec      = vcat(gasspec,solidspec)

# Assign numbers to each species based on their order in the merged list
for (i, specie) in enumerate(spec)
    symbol = Symbol(specie * "ₙ")
    @eval const $(symbol) = $i
end

f(α,β) = α*β

# species that has reactions "α", "β", "γ",  "θ", "ζ", "η", "ϕ" 
# species Ids αₙ, βₙ, γₙ, θₙ and so on so forth 
# reaction of α is a function f(α,β) = α*β and the function arguments IDs are [αₙ, βₙ]
# reaction of β is a function f(θ,β) = θ*β and the function arguments IDs are [θₙ, βₙ]
# reaction of γ is a function f(γ,ϕ) = γ*ϕ and the function arguments IDs are [γₙ, ϕₙ]
