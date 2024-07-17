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
const Eb  = 135e3J/mol # Activation enegy of pyrolysis by kissinger method  Unit: J.mole -1
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
const νc   = 0.55    # char yield coefficient                        Unit: -
const νO2  = 1.7      # oxygen yeild coefficient                      Unit: -                       

# Stoichiometric Parameters
const Mg  = 0.02897kg/mol             # Molar Weight of Air         Unit: Kg mol -1
const Yb0  = 1                   # Initial bitumen concentration    Unit: -
const Yc0  = 0                   # Initial char concentration          Unit: -
#const Pg   = atm           # Gas Pressure                      Unit: Pa
P = atm + 50
# Initial Temperatures and density
const Ta = 298K        # Ambient temperature;                               Unit: K
const ρᵧ₀   = (Mg .* P) / (R .* Ta)# Initial Density of Air (Gas)Unit: kg m -3
const YO20 = 0.204   
const Rgas = R/Mg 

th = 4865; # heater off time 
tg = 4532; # Air on time 

grid = simplexgrid(0:h:Lx)
gridplot(grid;Plotter=Plotter,resolution = (600, 250), legend =:rt)
using LaTeXStrings
species = ["Yᵦ", "Yᵪ", "ρᵧ", "Yₒ", "Tₛ", "Tᵧ"]
SpecInd = []
for (i, specie) in enumerate(species)
    symbol = Symbol(specie * "ₙ")
    @eval const $(symbol) = $i
	push!(SpecInd,i)
end

#=
mutable struct StorageStruct
    SpeciesId::Int
    Func::Function
    Const
    ArgsIds::Vector{Int}
end
=#

#=
StoreList = [StorageStruct(Yᵦₙ, ()->[], 1, []),StorageStruct(Yᵪₙ, ()->[], 1, []), StorageStruct(Pₙ, (Tᵧ)->ϕg/Rgas/Tᵧ, 1, [Tᵧₙ]),
StorageStruct(Yₒₙ, (Tᵧ)->ϕg/Rgas/Tᵧ, 1, [Tᵧₙ])
] =#



["Yᵦ", "Yᵪ", "P", "Yₒ", "Tₛ", "Tᵧ"]
function storagesmold!(f,u,node,data)
    f[Yᵦₙ] = u[Yᵦₙ]
	f[Yᵪₙ] = u[Yᵪₙ]
	f[ρᵧₙ] = ϕg*u[ρᵧₙ]
	f[Yₒₙ] = ϕg*u[Yₒₙ]
	f[Tₛₙ] = ρCpeff*u[Tₛₙ]
	f[Tᵧₙ] = ϕg*Cpg*u[ρᵧₙ]*u[Tᵧₙ]
end


function fluxsmold!(f,u,edge,data)

	# Flux for ρg
    nk = edge.node[1]
    nl = edge.node[2]
    Δx = edge.coord[nl]-edge.coord[nk]
    # arithmetic mean 
   
    # Gas velocity 
    uᵧ         = (Permiabilit/μ)*R*(u[Tᵧₙ,1]*u[ρᵧₙ,1]-u[Tᵧₙ,2]*u[ρᵧₙ,2])
	# convective upwind flux in ρᵧ: which is gas density 
	f[ρᵧₙ]     = uᵧ>0 ? uᵧ*u[ρᵧₙ,1] : uᵧ*u[ρᵧₙ,2] 
	# Oxygen species transport convective upwind flux
	ρₘ 		   = (u[ρᵧₙ,1] + u[ρᵧₙ,2]) # Arithmetic mean of density for diffusion and velocity calculations.. 
	convUpwind = ρₘ*uᵧ>0 ? ρₘ*uᵧ*u[Yₒₙ,1] : ρₘ*uᵧ*u[Yₒₙ,2] 
	# total species flux
	f[Yₒₙ] 		= convUpwind+ϕg*Dg*ρₘ*(u[Yₒₙ,1]-u[Yₒₙ,2])
	# diffusive flux for solid temperature 
	f[Tₛₙ] 		= keff*(u[Tₛₙ,1]-u[Tₛₙ,2])
	# Flux for gas temperature 
	uvelTgas 	= ρₘ*Cₚᵧ*uᵧ
	convUpwindT = uvelTgas>0 ? uvelTgas*u[Tᵧ,1] : uvelTgas*u[Tᵧ,2] 
	f[Tᵧ] 		= convUpwindT+ϕg*kᵧ*(u[Tᵧ,1]-u[Tᵧ,2])
end

# ["Yᵦₙ", "Yᵪₙ", "ρᵧₙ", "Yₒₙ", "Tₛₙ", "Tᵧₙ"]

function reactionsmold!(f,u,node,data)	
	Rb = Ab*exp(-Eb/(R*u[Tₛₙ]))*u[Yᵦₙ]
	Rc = Ac*exp(-Eb/(R*u[Tₛₙ]))*u[Yᵪₙ]*u[Yₒₙ]
	f[Yᵦₙ] = -Rb
	f[Yᵪₙ] = νc*Rb-Rc

	f[ρᵧₙ] = (ϕᵦ*ρᵦ)*((1-νc)*Rb+(1-νO2)*Rc)

	f[Yₒₙ] = -(ϕᵦ*ρᵦ)*νO2*Rc

	external = U*(AV_c)*(u[Tₛₙ]-Ta)
	internal = hsg*(Av_s)*(u[Tᵧₙ]-u[Tₛₙ])
	Q        = (ϕᵦ*ρᵦ)*(Hc*Rc+Hb*Rb)

	f[Tₛₙ] = -external + internal - Q
	f[Tᵧₙ] = -internal
end

function outflow!(f, u, node,data)
    if node.region == 1
		out = ugas * u[ρᵧ];
        f[ρᵧ] = ramp(node.time, dt = (0.0,tg), du = (0.0,out))
    end
end

function bconditionsmold!(f,u,bnode,data)

	# boundary condition for the density 
	#value=P/Rgas/T0
	#v=ramp(bnode.time,dt=(0,tg),du=(value, 0.0))
	#boundary_neumann!(f,u,bnode,species=ρᵧ,region=1,v)
	boundary_dirichlet!(f,u,bnode,species=ρᵧ,region=2,value=P/Rgas/T0)


	# Boundary condition for YO2, Sp = 4
	boundary_dirichlet!(f,u,bnode,species = Yₒₙ,region=1,value=YO20)
	boundary_robin!(f, u,bnode,species = Yₒₙ, region = 2, factor = hm, value = hm*YO20 )
	#boundary_neumann!(f,u,bnode,species=5,region=2,value=0) # Neumann for Ts at x = L
	
	# Boundary condition for Ts, Sp = 5
	v=ramp(bnode.time,dt=(0,th),du=(25000/keff, 0.0))
	boundary_neumann!(f,u,bnode,species = Tₛₙ,region=1,value=v)
	boundary_neumann!(f,u,bnode,species = Tₛₙ,region=2,value=0) # Neumann for Ts at x = L
	
	# Boundary condition for Tg, Sp = 6
	boundary_dirichlet!(f,u,bnode,species = Tᵧₙ,region=1,value=298)
	boundary_neumann!(f,u,bnode,species = Tᵧₙ,region=2,value=0) # Neumann for Tg at x = L
end


SmoldSystem  =   VoronoiFVM.System(grid;
    data       = data,
    flux       = fluxsmold!,
    storage    = storagesmold!,
    reaction   = reactionsmold!,
	bcondition = bconditionsmold!,
    unknown_storage = :dense,
    assembly = :cellwise)
# activate the species in the fuel electrode: Region: 1

for i in SpecInd
    enable_species!(SmoldSystem, SpecInd[i], [1])	
end



inivalSmold=unknowns(SmoldSystem)

inivalSmold[Tₛₙ,:] .= Ta
inivalSmold[Tᵧₙ,:] .= Ta
inivalSmold[Yₒₙ,:] .= YO20
inivalSmold[ρᵧₙ,:] .= ρᵧ₀
inivalSmold[Yᵦₙ,:] .= 1.0
inivalSmold[Yᵪₙ,:] .= 0.0



smold_ode_problem = ODEProblem(SmoldSystem,inivalSmold,(0,final_t))

Sol = solve(
    smold_ode_problem,
		ImplicitEuler(),
		adaptive=true,
        reltol=1.0e-4,
		abstol=1.0e-4,
		initializealg=NoInit()
	)	
smoldsol=reshape(Sol,SmoldSystem)

	