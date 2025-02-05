using LinearAlgebra, Printf, LaTeXStrings
using ExtendableGrids, VoronoiFVM
using ExtendableSparse, LinearSolve, ForwardDiff
using GridVisualize, GLMakie
using Plots, OrdinaryDiffEq, DataFrames, LessUnitful
using SimplexGridFactory
using Triangulate
using TetGen
using CairoMakie
using ColorSchemes

# ---- Physical Parameters ----
const ρ_c = 1500.0       # kg/m^3 (Coal density)
const ρₐ = 1.1          # kg/m^3 (Air density)
const C_pc = 1000.0      # J/kg-K (Specific heat of coal)
const C_pa = 1000.0      # J/kg-K (Specific heat of air)
const Dₐ = 2e-5         # m^2/s (Diffusion coefficient of oxygen in air)
const ε = 0.2           # Porosity
const Eₐ = 7e4          # J/mol (Activation energy)
const Bₐ = 4.2e-10      # 1/s (Pre-exponential factor)
const R = 8.314         # J/mol-K (Universal gas constant)
#const V = 0.6e-5      # m/s (Velocity in x-direction)
const λₛ = 0.12         # W/m-K (Thermal conductivity)
const ΔH = 3e5          # J/mol-O₂ (Heat of reaction)
const MW_O2 = 32e-3     # kg/mol (Molecular weight of O₂)
const T₀ = 293.0        # K (Initial temperature)
const c1₀ = 0.21        # Initial oxygen concentration
const n = 1.0           # Reaction order
const f_c2 = 1.0        # Weathering function
const Days = 350        # Simulation time (days)

# ---- Grid Parameters (2D) ----
const Lx = 10.0   # Length of coal pile in x-direction (m)
const Ly = 5.0    # Length in y-direction (m)
const Nx = 40     # Number of grid points in x
const Ny = 20     # Number of grid points in y
const Δt = 7200.0  # Time step (s)

# ---- Velocity (X and Y Components) ----
const V_x = 0  # m/s (Velocity in x-direction)
const V_y = -0.6e-5  # m/s (Smaller velocity in y-direction)
const Vr = (V_x, V_y)  # Velocity as a vector
# ---- Create 2D Grid ----
#grid = simplexgrid(0:Lx/Nx:Lx, 0:Ly/Ny:Ly)
#gridplot(grid; Plotter=GLMakie, resolution=(600, 400), legend=:rt)

builder_coal = let
    b = SimplexGridBuilder(; Generator = Triangulate)

    # Define points for a hill with a flat top
    p1 = point!(b, 0.0, 0.0)   # Left bottom
    p2 = point!(b, 4.0, 10.0)  # Left peak
    p3 = point!(b, 6.0, 10.0)  # Right peak (flat top starts)
    p4 = point!(b, 10.0, 0.0)  # Right bottom

    # Define outer boundary with facets
    facetregion!(b, 2)
    facet!(b, p1, p2)
    #facetregion!(b, 2)
    facet!(b, p2, p3)
    #facetregion!(b, 4)
    facet!(b, p3, p4)
    facetregion!(b, 1)
    facet!(b, p4, p1)

    # Set mesh refinement - smaller value for better quality
    options!(b; maxvolume = 1/10)

    # Return the builder
    b
end

# Generate the grid from the builder
coal_grid = simplexgrid(builder_coal)

# Ensure visualization
vis = GridVisualizer(; resolution = (800, 800), Plotter = GLMakie)
gridplot!(vis, coal_grid; linewidth = 0.5, alpha = 0.8)  # Increase transparency to see details
reveal(vis)
println("Coal stockpile 2D grid successfully generated and visualized!")
# ---- Define Variables ----
c₁ = 1; T = 2; c₂ = 3
SpecInd = [c₁, T, c₂]

# ---- Storage Function ----
function storage!(f, u, node, sim_data)
    f[c₁] = ε * ρₐ * u[c₁]                   # Oxygen storage
    f[T] = (1 - ε) * ρ_c * C_pc * u[T]       # Temperature storage
    f[c₂] = u[c₂]                            # Adsorbed oxygen storage
    return nothing
end

#= ---- Flux Function ----
function flux!(f, u, edge, sim_data)
    # Project values correctly across the edge
    vel_c1 = project(edge, ρₐ * V )  # Scale velocity by edge normal in 2D
    vel_T = project(edge, ρₐ * C_pa * V )

    # Compute Bernoulli flux factors for convection and diffusion
    Bplus = ε * Dₐ * ρₐ * fbernoulli(vel_c1 / (ε * Dₐ * ρₐ))
    Bminus = ε * Dₐ * ρₐ * fbernoulli(-vel_c1 / (ε * Dₐ * ρₐ))
    #Bplus = ε * Dₐ * ρₐ * fbernoulli(ρₐ *  V / (ε * Dₐ * ρₐ))
    #Bminus = ε * Dₐ * ρₐ * fbernoulli(-ρₐ *  V / (ε * Dₐ * ρₐ))

    # Correct flux computation with proper edge indexing
    f[c₁] = Bminus * u[c₁, 1] - u[c₁, 2]

    # Compute temperature flux similarly
    Bplus_T = λₛ * fbernoulli(vel_T / λₛ)
    Bminus_T = λₛ * fbernoulli(-vel_T / λₛ)
    #Bplus_T = λₛ * fbernoulli(ρₐ *  C_pa *  V / λₛ)
    #Bminus_T = λₛ * fbernoulli(-ρₐ *  C_pa *  V / λₛ)
    f[T] = Bminus_T * u[T, 1] - Bplus_T * u[T, 2]

    return nothing
end
=#

real_value(x) = x isa ForwardDiff.Dual ? ForwardDiff.value(x) : x

function flux!(f, u, edge, sim_data)
    # Define velocity components in x and y direction
    V_x = 0  # Velocity in x-direction (m/s)
    V_y = -0.6e-4  # Velocity in y-direction (m/s), less than x to mimic real-life scenario

    # Define velocity vector
    velocity_vector = [V_x, V_y]

    # Project velocity onto the edge (accounts for edge direction)
    V_proj = project(edge, velocity_vector)  

    # Extract values properly from dual numbers
    u_c1_edge = [real_value(u[c₁, i]) for i in 1:2]  # Extract for both edge points
    u_T_edge = [real_value(u[T, i]) for i in 1:2]

    # Compute Bernoulli flux factors for convection and diffusion
    Bplus = ε * Dₐ * ρₐ * fbernoulli(V_proj / (ε * Dₐ * ρₐ))
    Bminus = ε * Dₐ * ρₐ * fbernoulli(-V_proj / (ε * Dₐ * ρₐ))

    # Compute oxygen flux
    f[c₁] = Bminus * u_c1_edge[1] - Bplus * u_c1_edge[2]

    # Compute temperature flux
    Bplus_T = λₛ * fbernoulli(V_proj / λₛ)
    Bminus_T = λₛ * fbernoulli(-V_proj / λₛ)

    f[T] = Bminus_T * u_T_edge[1] - Bplus_T * u_T_edge[2]

    return nothing
end

# ---- Reaction Function ----
function reaction!(f, u, node, sim_data)
    r_ads = -Bₐ * u[c₁]^n * f_c2 * abs(exp(-Eₐ / R * (1/u[T] - 1/T₀)))
    f[c₁] = - (1 - ε) * ρ_c * r_ads
    f[T] = (1 - ε) * ρ_c * ΔH / MW_O2 * r_ads
    f[c₂] = r_ads
    return nothing
end

# ---- Boundary Conditions ----
function bcondition!(f, u, bnode, sim_data)
    boundary_dirichlet!(f, u, bnode, species=c₁, region=1, value=0.0)  # No oxygen at the bottom
    boundary_dirichlet!(f, u, bnode, species=T, region=1, value=T₀)  # Fixed temperature at the bottom

    boundary_dirichlet!(f, u, bnode, species=c₁, region=2, value=c1₀)  # Oxygen supply at the top
    boundary_dirichlet!(f, u, bnode, species=T, region=2, value=T₀)  # Fixed temperature at the top
end

# ---- Create System ----
sim_data = (
    Bₐ = Bₐ,
    C_pc = C_pc,
    C_pa = C_pa,
    Dₐ = Dₐ,
    ε = ε,
    Eₐ = Eₐ,
    R = R,
    Vr = Vr,
    λₛ = λₛ,
    ΔH = ΔH,
    MW_O2 = MW_O2,
    T₀ = T₀,
    c1₀ = c1₀,
    n = n,
    f_c2 = f_c2,
)

sys = VoronoiFVM.System(coal_grid; data=sim_data, 
                        flux=flux!, 
                        storage=storage!, 
                        reaction=reaction!, 
                        bcondition=bcondition!, 
                        unknown_storage=:sparse, 
                        assembly=:edgewise, species=SpecInd)

# ---- Initial Conditions ----
inival = unknowns(sys)
inival[c₁, :] .= c1₀  
inival[T, :] .= T₀    
inival[c₂, :] .= 0.0  

# ---- Time Evolution ----
tend = Days * 24 * 60 * 60
problem = ODEProblem(sys, inival, (0.0, tend))
sol = solve(problem, Rosenbrock23(); dt=Δt, reltol=1e-5, abstol=1e-5)
sol = reshape(sol, sys)

# Ensure visualization setup
vis_T = GridVisualizer(; size = (500, 600), Plotter=GLMakie)
vis_c1 = GridVisualizer(; size = (500, 600), Plotter=GLMakie)
vis_c2 = GridVisualizer(; size = (500, 600), Plotter=GLMakie)

# Plot Temperature Distribution
scalarplot!(vis_T, coal_grid, sol[end][T, :];
    title = "Final Temperature Distribution",
    colormap = :jet,  
    contour_lines = false,
    colorbar_label = "Temperature (K)",
    clear = false,
    xlabel="Distance", ylabel="Temperature (K)"
)

# Plot Oxygen Concentration
scalarplot!(vis_c1, coal_grid, sol[end][c₁, :];
    title = "Final Oxygen Concentration",
    colormap = :jet,
    contour_lines = false,
    colorbar_label = "Oxygen Concentration (%)",
    clear = false,
    xlabel="Distance", ylabel="Oxygen Concentration (%)"
)

# Plot Adsorbed Oxygen
scalarplot!(vis_c2, coal_grid, sol[end][c₂, :];
    title = "Final Adsorbed Oxygen",
    colormap = :jet,
    contour_lines = false,
    colorbar_label = "Oxygen adsorbed in coal (kg O₂/kg coal)",
    clear = false,
    xlabel="Distance", ylabel="Oxygen adsorbed in coal (kg O₂/kg coal)"
)

# Display plots
reveal(vis_T)
reveal(vis_c1)
reveal(vis_c2)

#= Ensure visualization setup
vis_T = GridVisualizer(; size = (800, 800), Plotter=GLMakie)
vis_c1 = GridVisualizer(; size = (800, 800), Plotter=GLMakie)
vis_c2 = GridVisualizer(; size = (800, 800), Plotter=GLMakie)

# Plot Temperature Distribution
scalarplot!(vis_T, coal_grid, sol[end][T, :];
title = "Final Temperature Distibution",
colormap = :jet,  
contour_lines = :false,
colorbar_title = "Temperature (°C)"  
)

# Plot Oxygen Concentration
scalarplot!(vis_c1, coal_grid, sol[end][c₁, :]; title="Final Oxygen Concentration", colormap=:jet, contour_lines = false)

# Plot Adsorbed Oxygen
scalarplot!(vis_c2, coal_grid, sol[end][c₂, :]; title="Final Adsorbed Oxygen", colormap=:jet, contour_lines = false)

# Display plots
reveal(vis_T)
reveal(vis_c1)
reveal(vis_c2)
=#