using firesolid

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
using GLMakie
using Colors

#using Plots
Plotter = GLMakie
lib_dir = "lib/"
tr_file = joinpath(lib_dir, "transport.dat")
therm_file = joinpath(lib_dir, "therm.dat")
@unitfactors mol dm eV μA μF cm μm bar atm R m s A Pa mm kg J K W
@phconstants N_A e;

const Lx  = 0.45m;  # Length of column above the heater (Bitumen+sand)    Unit: m
const Ls  = 0.1m;   # Length of column below the heater (clean sand)      Unit: m
const h  = 0.01m;   # grid size 
const dt      = 20s;          # Time step
const final_t = 20000s;       # Total Time (sec)


# Geometric Parameters
const rc  = 0.08m    # Radius of the column   Unit: m
const dp  = 0.00088m  # Particle diameter size  Unit: m
const N   = 11       # Number of thermocouples    
const Ni  = 0.12 - Ls  # Location of first thermocouple(0.12m absolute)   
const ϕₛ  = 0.055      # Sand porosity                                    
const ϕᵦ  = ϕₛ         # Bitumen porosity                                 
const ϕg  = 0.315      # Gas filled porosity                              
const ϕₜ  = ϕₛ + ϕg  # Total porosity                                   
const S_A_c = 2*π*rc*(Lx-Ls) # Surface area of column above heater 
const V_c = π*rc^2*(Lx-Ls) # Volume of column above heater (Bitumen+sand)Unit: m3
const AV_c = S_A_c/ V_c    # Surface area to volume ratio for a cylinder 

const S_A_s = 4 * π * (dp / 2)^2           # Surface area of sphere (Sand 
const V_s = (4 / 3) * π * (dp / 2)^3     # Volume of sphere (Sand Particle)      
const AV_s = (1-ϕₜ)*S_A_s/V_s # Surface area to volume ratio for Thermal Parameters


const ρₛ = 2650kg/m^3               # Density of sand particle      Unit: kg.m-3         
const ρᵦ = 1030kg/m^3               # Density of bitumen     Unit: kg.m-3  
const ρ  = 1726kg/m^3               # Bulk density of combined sand and bitumen        
const Cpb = 921J/kg/K                 # Bitumen specific heat capacity  Unit: J kg −1 K −1
const kb  = 0.15W/m/K                # Thermal conductivity of bitumen Unit: W m −1 K −1
const U   = 13W/m^2/K                   # Global heat loss coefficient Unit: W m −2 K −1
const σ   = 5.67e-8W/m^2/K^4          # Stefan boltzman constant Unit: W m −2 K −4
const Ac  = 10^4.9          # Arhenius coefficient for char oxidation(log(Ac)= 4.9)
const Ab  = 10^7.5 # Arhenius coefficient for bitumen pyrolysis       (log(Ab)= 7.5)
const Dg  = 4.53e-5m^2/s    # Diffusion coefficient for gas                Unit: m 2 s −1
const Pr  = 0.72                 # Prandtl Number
const Re  = 11                   # Reynolds number for typical smouldering
const Ep  = 135e3J/mol # Activation enegy of pyrolysis by kissinger method  Unit: J.mole -1
const Ec  = 90e3J/mol # Activation enegy of oxidation of char     Unit: J.mole -1
const Hb  = 1.62e6J/kg # Heat of reaction for pyrolysis (Endothermic)   Unit: J kg −1
const Hc  = -38.73e6J/kg # Heat of reaction for oxidation (Exothermic)      Unit: J kg −1
const hm  = 100kg/m^2/s # Mass transfer coefficient of O2       Unit: kg m −2 s −1             
#const Ru  = 8.3145               # Universal gas constant   Unit: J K-1 mol-1
const μ = 1.80e-5kg*s/m^3  # Viscosity of air at room temperature Unit: kg s m -3  (Pa.s)
const κ = 1e-9m^2  # Permiability of the porous media                   Unit: m2

const vf = 6.933e-5m/s # Smouldering front end velocity                Unit: m s -1                                
const ug_all  = [0.058, 0.025, 0.05, 0.08]m/s   # Darcy air flux Unit: m s -1
const ug = ug_all[1]    # Darcy air flux considered in this work        Unit: m s -1
const vc   = 0.55    # char yield coefficient                        Unit: -
const vO2  = 1.7      # oxygen yeild coefficient                      Unit: -                       

# Stoichiometric Parameters
const Mg  = 0.02897kg/mol             # Molar Weight of Air         Unit: Kg mol -1
const Yb0  = 1                   # Initial bitumen concentration    Unit: -
const Yc0  = 0                   # Initial char concentration          Unit: -
#const Pg   = atm           # Gas Pressure                      Unit: Pa

# Initial Temperatures and density
const Ta = 298K        # Ambient temperature;                               Unit: K
const pho_g0   = (Mg .* atm) / (R .* Ta)# Initial Density of Air (Gas)Unit: kg m -3
const YO20 = 0.204   
const Rgas = R/Mg 

grid = simplexgrid(0:h:Lx)
gridplot(grid;Plotter=Plotter,resolution = (600, 250), legend =:rt)
using LaTeXStrings
species = ["Yᵦ", "Yᵪ", "ρᵧ", "Yₒ", "Tₛ", "Tᵧ"]

for (i, specie) in enumerate(species)
    symbol = Symbol(specie * "ₙ")
    @eval const $(symbol) = $i
end

#=
mutable struct StorageStruct
    SpeciesId::Int
    Func::Function
    Const
    ArgsIds::Vector{Int}
end
=#

StoreList = [StorageStruct(Yᵦₙ, ()->[], 1, []),StorageStruct(Yᵪₙ, ()->[], 1, []), StorageStruct(Pₙ, (Tᵧ)->ϕg/Rgas/Tᵧ, 1, [Tᵧₙ]),
StorageStruct(Yₒₙ, (Tᵧ)->ϕg/Rgas/Tᵧ, 1, [Tᵧₙ])



]



["Yᵦ", "Yᵪ", "P", "Yₒ", "Tₛ", "Tᵧ"]
function storage(f,u,node,data)
    f[Yᵦₙ] = u[Yᵦₙ]
	f[Yᵪₙ] = u[Yᵪₙ]
	f[ρᵧₙ] = ϕg*u[ρᵧₙ]
	f[Yₒₙ] = ϕg*u[Yₒₙ]
	f[Tₛₙ] = ρCpeff*u[Tₛₙ]
	f[Tᵧₙ] = ϕg*Cpg*u[ρᵧₙ]*u[Tᵧₙ]
end


function flux(f,u,edge,data)

	# Flux for ρg
    nk = edge.node[1]
    nl = edge.node[2]
    Δx = edge.coord[nl]-edge.coord[nk]
    # arithmetic mean 
   
    
    uᵧ = - (Permiabilit/μ)*R*(u[Tᵧₙ,1]*u[ρᵧₙ,1]-u[Tᵧₙ,2]*u[ρᵧₙ,2])/Δx

	f[ρᵧₙ]=(Rgas*κ/μ)*umean*(u[Tᵧₙ,1]*u[ρᵧₙ,1]-u[Tᵧₙ,2]*u[ρᵧₙ,2]) ## check this for u[3] term
	
	# Oxygen species transport has two fluxes 
	flux1 = u[Yₒₙ]*(Rgas*κ/μ)*(u[ρᵧₙ,1]*u[Tᵧₙ,1]-u[ρᵧₙ,2]*u[Tᵧₙ,2])
	flux2 = ϕg*Dg*(u[4,1]-u[4,2])
	f[4] = flux1+flux2

	
	f[5] = keff*(u[5,1]-u[5,2])
	
	fluxT1 = u[3]*Cpg*u[6]*(R*κ/μ)*(u[3,1]*u[6,1]-u[3,2]*u[6,2])
	
	f[6] =fluxT1+ϕg*kg*(u[6,1]-u[6,2])
end