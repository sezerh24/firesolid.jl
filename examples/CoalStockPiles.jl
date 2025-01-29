using LinearAlgebra, Printf, LaTeXStrings
using ExtendableGrids, VoronoiFVM
using ExtendableSparse, LinearSolve
using GridVisualize
#using CairoMakie
using GLMakie
using Plots
using OrdinaryDiffEq
using DataFrames
using LessUnitful

# Physical Parameters
const ρ_c = 1500.0         # kg/m^3
const ρₐ = 1.1             # kg/m^3
const C_pc = 1000.0        # J/kg-K
const C_pa = 1000.0        # J/kg-K
const Dₐ = 2e-5            # m^2/s
const ε = 0.2              # Porosity
const Eₐ = 7e4             # J/mol
const Bₐ = 4.8e-10         # 1/s
const R = 8.314            # J/mol-K
const V = 0.6e-5           # m/s
const λₛ = 0.12             # W/m-K
const ΔH = 3e5             # J/mol-O2
const MW_O2 = 32e-3        # kg/mol
const T₀ = 293.0           # Initial temperature (K)
const c1₀ = 0.21           # Initial oxygen concentration
const n = 1.0              # Reaction order
const f_c2 = 1.0           # Weathering function
const Days = 472           # Simulation time (days)

# Additional Parameters
const L = 10.0             # Length of the coal pile (m)
const Nx = 60              # Number of grid points
const Δt = 900.0           # Time step (s)


# Create DataFrame from the parameters
tabledata = DataFrame(
    Property = [
        "Pre-exponential factor, Bₐ (1/s)",
        "Specific heat (coal), C_pc (J/kg-K)",
        "Specific heat (air), C_pa (J/kg-K)",
        "Diffusion coefficient (air), Dₐ (m^2/s)",
        "Porosity, ε",
        "Activation energy, Eₐ (J/mol)",
        "Gas constant, R (J/mol-K)",
        "Velocity, V (m/s)",
        "Thermal conductivity (solid), λₛ (W/m-K)",
        "Heat of reaction, ΔH (J/mol-O₂)",
        "Molecular weight of O₂, MW_O2 (kg/mol)",
        "Initial temperature, T₀ (K)",
        "Initial oxygen concentration, c1₀",
        "Reaction order, n",
        "Weathering function, f_c2"
    ],
    Value = [
        Bₐ,       # Pre-exponential factor
        C_pc,     # Specific heat of coal
        C_pa,     # Specific heat of air
        Dₐ,       # Diffusion coefficient
        ε,        # Porosity
        Eₐ,       # Activation energy
        R,        # Gas constant
        V,        # Velocity
        λₛ,       # Thermal conductivity
        ΔH,       # Heat of reaction
        MW_O2,    # Molecular weight of oxygen
        T₀,       # Initial temperature
        c1₀,      # Initial oxygen concentration
        n,        # Reaction order
        f_c2      # Weathering function
    ]
)              

c₁ = 1;
T = 2; 
c₂ = 3; 
SpecInd = [c₁,T,c₂]

# Create Grid
grid = simplexgrid(0:L/Nx:L)
gridplot(grid;Plotter=GLMakie,resolution = (600, 250), legend =:rt)

function storage!(f, u, node, data)
    f[1] =  ε * ρₐ * u[c₁]                     # Storage for oxygen concentration
    f[2] = (1 -  ε) * ρ_c *  C_pc * u[T]   # Storage for temperature
    f[3] = u[c₂]                                   # Storage for adsorbed oxygen
    return nothing
end

function flux!(f, u, edge, data)
    # Convective and diffusive flux for c1
    vel_c1 = project(edge, ρₐ *  V)
    Bplus =  ε *  Dₐ * ρₐ * fbernoulli(vel_c1 / ( ε *  Dₐ * ρₐ))
    Bminus =  ε *  Dₐ * ρₐ * fbernoulli(-vel_c1 / ( ε *  Dₐ * ρₐ))
    f[1] = Bminus * u[c₁, 1] - Bplus * u[c₁, 2]

    # Convective and diffusive flux for T
    vel_T = project(edge, ρₐ *  C_pa *  V)
    Bplus =  λₛ * fbernoulli(vel_T /  λₛ)
    Bminus =  λₛ * fbernoulli(-vel_T /  λₛ)
    f[2] = Bminus * u[T, 1] - Bplus * u[T, 2]

    return nothing
end

function reaction!(f, u, node, data)
    # Reaction rate for adsorbed oxygen
    r_ads =  Bₐ * u[c₁]^ n *  f_c2 * exp(- Eₐ /  R * (1 / u[T] - 1 /  T₀))

    # Source term for oxygen concentration
    f[1] = - (1 -  ε) * ρ_c * r_ads

    # Source term for temperature
    f[2] = (1 -  ε) * ρ_c *  ΔH /  MW_O2 * r_ads

    # Source term for adsorbed oxygen
    f[3] = r_ads
    return nothing
end

function bcondition!(f, u, bnode, data)
    boundary_dirichlet!(f, u, bnode, species = c₁, region = 1, value =  c1₀)
    boundary_dirichlet!(f, u, bnode, species = c₁, region = 2, value = 0.0)
    boundary_dirichlet!(f, u, bnode, species = T, region = 1, value =  T₀)
    boundary_dirichlet!(f, u, bnode, species = T, region = 2, value =  T₀)
    return nothing
end

data = (
    Bₐ = Bₐ,
    C_pc = C_pc,
    C_pa = C_pa,
    Dₐ = Dₐ,
    ε = ε,
    Eₐ = Eₐ,
    R = R,
    V = V,
    λₛ = λₛ,
    ΔH = ΔH,
    MW_O2 = MW_O2,
    T₀ = T₀,
    c1₀ = c1₀,
    n = n,
    f_c2 = f_c2,
)

sys = VoronoiFVM.System(grid;data = data,
flux       = flux!,
storage    = storage!,
reaction   = reaction!,
bcondition = bcondition!,
unknown_storage = :sparse,
assembly = :edgewise, species = SpecInd)


inival = unknowns(sys)
inival[c₁, :] .=  c1₀  # Initial oxygen concentration
inival[T, :] .=  T₀    # Initial temperature
inival[c₂, :] .= 0.0       # Initial adsorbed oxygen


# Simulation Parameters
tend = Days * 24 * 60 * 60  # Total simulation time (s)
problem = ODEProblem(sys, inival, (0.0, tend))

# Solver
#sol = solve(problem, Rosenbrock23(); dt = Δt, reltol = 1e-5, abstol = 1e-5)
sol=solve(problem,ImplicitEuler(),  

                                   force_dtmin=true,
                                   adaptive=true,
                                   reltol=1.0e-3,
                                   abstol=1.0e-3,
                                   initializealg=NoInit(),
                                  dtmin=1e-6,
                                  # force_dtmin = true,
                                   #progress=true,
                                   #progress_steps=1,

                                   )

sol=reshape(sol,sys)

# Visualization setup
vis = GridVisualizer(Plotter = GLMakie, resolution = (1200, 400), layout = (1, 3), legend = :rb)

# Final time visualization
scalarplot!(vis[1, 1], grid, sol[end][c₁, :], label = "Oxygen Concentration (c₁)")
scalarplot!(vis[1, 2], grid, sol[end][T, :], label = "Temperature (T)")
scalarplot!(vis[1, 3], grid, sol[end][c₂, :], label = "Adsorbed Oxygen (c₂)")
reveal(vis)

# Temporal Evolution
time_points = [139, 278, 471] .* 24 .* 60 .* 60  # Days in seconds

for t in time_points
    # Plot oxygen concentration
    scalarplot!(vis[1, 1], grid, sol(t)[c₁, :], 
                label = "c₁ at day $(t / (24 * 60 * 60))", clear = false)

    # Plot temperature
    scalarplot!(vis[1, 2], grid, sol(t)[T, :], 
                label = "T at day $(t / (24 * 60 * 60))", clear = false)

    # Plot adsorbed oxygen
    scalarplot!(vis[1, 3], grid, sol(t)[c₂, :], 
                label = "c₂ at day $(t / (24 * 60 * 60))", clear = false)
end
