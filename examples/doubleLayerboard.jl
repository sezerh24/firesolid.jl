# Step 2: Define the species and their IDs
gasspec   = ["α", "β", "γ", "θ", "ζ", "η", "ϕ"]
solidspec = ["ξ", "μ"]
species_list = vcat(gasspec, solidspec)

for (i, specie) in enumerate(species_list)
    symbol = Symbol(specie * "ₙ")
    @eval const $(symbol) = $i
end

# Step 3: Define the reactions
reactions = [
    Reaction(αₙ, (α, β) -> α * β, [αₙ, βₙ]),
    Reaction(βₙ, (θ, β) -> θ * β, [θₙ, βₙ]),
    Reaction(γₙ, (γ, ϕ) -> γ * ϕ, [γₙ, ϕₙ])
]

storages = StorageStruct(αₙ, (α,β)->α*β, [], [αₙ,βₙ])

A = rand(Float64,(4,25))
storagesF = [
    StorageF(αₙ, (α, β) -> α * β, [αₙ, βₙ])
]

storagesC = [
    StorageC(γₙ, A[1,:])
]

data = (reactions = reactions, StC = storagesC, StF = storagesF)

