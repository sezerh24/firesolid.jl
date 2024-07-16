mutable struct Storage{T}
    species_id::Int
    data::T
end

# Define different types of storage
struct StorageFunc
    func::Function
    args_ids::Vector{Int}
end

struct StorageConst
    const_data::Vector{Float64}
end

# Example functions using multiple dispatch
function process_storage(s::Storage{StorageFunc})
    println("Processing StorageFunc with species ID $(s.species_id)")
    # Access fields of StorageFunc using s.data.func, s.data.args_ids, etc.
end

function process_storage(s::Storage{StorageConst})
    println("Processing StorageConst with species ID $(s.species_id)")
    # Access fields of StorageConst using s.data.const_data, etc.
end

# Creating instances of the structs
sf = Storage{StorageFunc}(1, StorageFunc(x -> 2x, [1, 2]))
sc = Storage{StorageConst}(2, StorageConst([3.0, 4.0]))


