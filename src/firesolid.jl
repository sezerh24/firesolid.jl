module firesolid

# Additional packages needed. 
using DiffusionFlux
using ExtendableGrids
using ExtendableSparse
using GridVisualize
using IdealGas
using LessUnitful
using LinearAlgebra
using LinearSolve
using OrdinaryDiffEq
using RxnHelperUtils
using StaticArrays
using TransportProperties
using VoronoiFVM
using DataFrames

# Add structures and the VoronoiFVM physics functions 
include("StructSystem.jl")

#Export structures 
export StorageStruct, ReactionStruct, SourceStruct, FluxStruct

# Export physics functions
export storage!, reaction!, source!, Fluxfunc!

# Export diffusion functions
export MultiComponentFickeanDiff, DGM_Diffusion, MultiComponentDGM

end
