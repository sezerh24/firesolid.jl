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
#using GLMakie
using Colors
using MAT
using SimplexGridFactory
using Triangulate
using TetGen
small = 1e-20;
using Plots
using DataFrames
using PrettyTables
using LessUnitful
Plotter = Plots
@unitfactors mol dm eV μA μF cm μm bar atm R m s Pa mm A
@phconstants N_A e
F = N_A*e # Faraday constant 
ΔV = -0.45
using DataFrames
using PrettyTables

# Creating the DataFrame with the given data
sim_param = DataFrame(
    Description = [
        "Exc. current density", "Surface tension", "Phase-field interface thickness",
        "Barrier height", "Gradient energy coefficient", "Anisotropy strength",
        "Anisotropy mode", "Kinetic coefficient", "Site density electrode",
        "Bulk Li-ion concentration", "Conductivity electrode", "Conductivity electrolyte",
        "Diffusivity electrode", "Diffusivity electrolyte"
    ],
    Symbol = [
        "i₀", "γ", "δ_{PF}", "W", "κ₀", "δ_{aniso}",
        "ω", "Lη", "C_{ms}",
        "C₀", "σₛ", "σₗ",
        "Dₛ", "Dₗ"
    ],
    RealValue = [
        "30 [A/m²]", "0.556 [J/m²]", "1.5 × 10⁻⁶ [m]",
        "4.45 × 10⁶ [J/m³]", "1.25 × 10⁻⁶ [J/m]", "0.044",
        "4", "1.81 × 10⁻³ [1/s]", "7.64 × 10⁴ [mol/m³]",
        "10³ [mol/m³]", "10⁷ [S/m]", "1.19 [S/m]",
        "7.5 × 10⁻¹³ [m²/s]", "3.197 × 10⁻¹⁰ [m²/s]"
    ],
    Normalized = [
        30, 0.22, 1.5,
        1.78, 0.5, 0.044,
        4, 1.81e-3, 76.4,
        1, 1.0e7, 1.19,
        0.75, 319.7
    ],
    Source = [
        "[123]", "[124, 125]", "[95]",
        "Computed", "Computed", "[97, 125]",
        "[97, 98]", "Computed", "[66]",
        "Computed", "[64]", "[126]",
        "[64]", "[126]"
    ]
)

# Display the DataFrame as a pretty table
#pretty_table(sim_param, header = ["Description", "Symbol", "Real Value", "Normalized", "Source"])

# assign SpecIDs 
ξₙ = 1;
ζₙ = 2;
ϕₙ = 3;
SpecInd = [ξₙ, ζₙ, ϕₙ]
# Simulation parameters 
i₀ = sim_param.Normalized[1]  # Exchange current density
γ = sim_param.Normalized[2]  # surface tension 
δpf = sim_param.Normalized[3] # Phase field interface thickness 
W = sim_param.Normalized[4]  # Barrier height
κ₀ = sim_param.Normalized[5] # Gradient energy coefficient 
δaniso = sim_param.Normalized[6] # Anisotropy strength 
ω   = sim_param.Normalized[7] #Anisotropy mode 
Lₙ = sim_param.Normalized[8] # Kinetic coefficient 
Cₘₛ = sim_param.Normalized[9] # Site density electrode 
C₀ = sim_param.Normalized[10] # Bulk Li-ion concentration 
σₛ = sim_param.Normalized[11] # Conductivity electrode 
σₗ = sim_param.Normalized[12] # Conductivity electrolyte 
Dₛ = sim_param.Normalized[13] # Diffusivity electrode 
Dₗ = sim_param.Normalized[14] # Diffusivity electrolyte 
α = 0.5; # charge transfer coefficient
ne = 2; # number of the electron transfer in the electrochemical reaction 
T       = 273+25; # K: temperature
Lσ = 2000; # interfacial mobility


K = (κ₀*(1+δaniso))
# Define the functions

#using Plots
h(ξ) = exp(20 * (ξ - 0.5)) / (1 + exp(20 * (ξ - 0.5))) 
g(ξ) = W * ξ^2 * (1 - ξ)^2
# Define the derivatives
dh_dξ(ξ) = 20 * h(ξ) * (1 - h(ξ)) 
dg_dξ(ξ) = W * 2 * ξ * (1 - ξ) * (1 - 2 * ξ)
#=
ξval = sort(rand(1000))
plot(ξval, h.(ξval))
plot!(ξval, dh_dξ.(ξval))

plot(ξ, g.(ξ))
plot!(ξ, dg_dξ.(ξ))
=#

#=
vis1=GridVisualizer(Plotter=Plotter;
resolution=(1200,300),
layout=(1,1), legend = :lt
)
=#

#scalarplot!(ξ,h.(ξ)')

Lx = 50
Ly = 20
Nx = 100
Ny = 40
X = 0:Lx/(Nx-1):Lx
Y = 0:Ly/(Ny-1):Ly
grid = simplexgrid(X,Y)
#gridplot(grid;Plotter=Plotter,resolution = (900, 450), legend =:rt)


function reactionLi!(f,u,node)  
   # f[ξₙ] = -Lₙ*dh_dξ(u[ξₙ])*((exp((1-α)*ne*F*u[ϕₙ]/R/T)) - u[ζₙ]*(exp(-α*ne*F*u[ϕₙ]/R/T))) - Lσ*dg_dξ(u[ξₙ]) 
   
    
    ξₙ_val = u[ξₙ]<=0 ? max(u[ξₙ],0.0) :  min(u[ξₙ],1.0)
    ζₙ_val  = u[ζₙ]<=0 ? max(u[ζₙ],0.0) :  min(u[ζₙ],1.0)
    ϕ_val = u[ϕₙ]>=0 ? min(0,u[ϕₙ]) : max(u[ϕₙ], ΔV)
    
#=
     arg1      = (exp((1-α)*ne*F*u[ϕₙ]/R/T)) 
     arg2      = (u[ζₙ])*abs(exp(-α*ne*F*u[ϕₙ]/R/T)) 
     Ln   = Lₙ*dh_dξ(abs(u[ξₙ])) 
     Lsig = Lσ*dg_dξ(abs(u[ξₙ]))
=#

     arg1      = exp((1-α)*ne*F*ϕ_val/R/T) 
     arg2      = ζₙ_val*exp(-α*ne*F*ϕ_val/R/T)
     Ln    = Lₙ*dh_dξ(ξₙ_val) 
     Lsig  =  Lσ*dg_dξ(ξₙ_val) 

    f[ξₙ] = -(Lsig+Ln*(arg1- arg2))
    f[ζₙ] = (Lsig+Ln*(arg1- arg2))*(Cₘₛ/C₀)
    f[ϕₙ] = (Lsig+Ln*(arg1- arg2))*ne*F*Cₘₛ
   #println("I am reaction")
#=
    if isnan(u[ξₙ].value) || isnan(u[ϕₙ].value) || isnan(u[ζₙ].value)
        println("One of the input values is NaN.")
    end
    =#
  #  f[ζₙ] = 0.0+ small
  #  f[ϕₙ] = 0.0+small
end

function storageLi!(f,u,node)

   ξₙ_val = u[ξₙ]<=0 ? max(u[ξₙ],0.0) :  min(u[ξₙ],1.0)
   ζₙ_val = u[ζₙ]<=0 ? max(u[ζₙ],0.0) :  min(u[ζₙ],1.0)
   ϕ_val = u[ϕₙ]>=0 ? min(0,u[ϕₙ]) : max(u[ϕₙ], ΔV)

    f[ξₙ] = ξₙ_val  
    f[ζₙ] = ζₙ_val 
    #f[ϕₙ] = ϕ_val

    
    #+ (Cₘₛ/C₀)*ξₙ_val 
    #f[ϕₙ] = u[ϕₙ]
    #f[ϕₙ] = ne*F*Cₘₛ*f[ξₙ] 
#=
    f[ξₙ] = u[ξₙ]  
    f[ζₙ] = u[ζₙ]
=#
    #=
    f[ξₙ] = u[ξₙ]  
    f[ζₙ] = u[ζₙ] +(Cₘₛ/C₀)*u[ξₙ]#ξₙ_val 
    f[ϕₙ] = ne*F*Cₘₛ*u[ξₙ]
    =#
   

     # Compute intermediate values and check for NaNs
     #ξₙ_val = f[ξₙ]
    # ζₙ_val = f[ζₙ]
     #ϕₙ_val = f[ϕₙ]
 
    #println("ξₙ: ", f[ξₙ], ", ζₙ: ", f[ζₙ], ", ϕₙ: ", f[ϕₙ])
    #println("I am storage")
end

function fluxLi!(f, u, edge)

    ξₘ = ((u[ξₙ,1] +u[ξₙ,2])/2.0)
    ζₘ = ((u[ζₙ,1] +u[ζₙ,2])/2.0)

  Deff = Dₛ*h(ξₘ)+Dₗ*(1-h(ξₘ)) 
  σeff = σₛ*h(ξₘ)+σₗ*(1-h(ξₘ)) 

  # calculate the diffusion fluxes 
  f[ξₙ] = Lσ*K*(u[ξₙ,1]-u[ξₙ,2]) 
  f[ζₙ] = Deff*((u[ζₙ,1]-u[ζₙ,2])) + (Deff*ne*F/(R*T))*ζₘ*(u[ϕₙ,1]-u[ϕₙ,2]) - (Cₘₛ/C₀)*Lσ*K*(u[ξₙ,1]-u[ξₙ,2])
  f[ϕₙ] = σeff*(u[ϕₙ,1]-u[ϕₙ,2]) - ne*F*Cₘₛ*Lσ*K*(u[ξₙ,1]-u[ξₙ,2])
   
  #=
  Vel = (Deff*ne*F/(R*T))*(u[ϕₙ,1]-u[ϕₙ,2]) 
  Bplus =  Deff*fbernoulli(Vel/Deff)
  Bminus =  Deff*fbernoulli(Vel/Deff)
  f[ζₙ] = Bminus * u[ζₙ,1] - Bplus * u[ζₙ,2] - (Cₘₛ/C₀)*Lσ*K*(u[ξₙ,1]-u[ξₙ,2])
=#
   
  #println(" I am flux ")
  #=
  if isnan(u[ξₙ].value) || isnan(u[ϕₙ].value) || isnan(u[ζₙ].value)
    println("One of the input values is NaN.")
  end
  =#

end



function bconditionsLi!(f,u,bnode)

    println("time = ", bnode.time)
	boundary_dirichlet!(f,u,bnode,species=ϕₙ,region=4,value=ΔV)
    boundary_dirichlet!(f,u,bnode,species=ϕₙ,region=2,value=0.0)
    boundary_dirichlet!(f,u,bnode,species=ξₙ,region=4,value=1.0)
    boundary_dirichlet!(f,u,bnode,species=ζₙ,region=4,value=0.0)
    boundary_dirichlet!(f,u,bnode,species=ζₙ,region=2,value=1.0) 
    #=
    if isnan(u[ξₙ]) || isnan(u[ϕₙ]) || isnan(u[ζₙ])
        println("One of the input values is NaN.")
    end
    =#
end


sysLi  =   VoronoiFVM.System(grid;
    flux       = fluxLi!,
    storage    = storageLi!,
    reaction   = reactionLi!,
	bcondition = bconditionsLi!,
    unknown_storage = :sparse,
    assembly = :cellwise)
# activate the species in the fuel electrode: Region: 1

for i in SpecInd
    enable_species!(sysLi, SpecInd[i], [1])	
end
inivalLi=float(unknowns(sysLi))


cξ = [0.0, Ly/2];
rx  = 4
ry  = 1
xξ(x,y)   = ((x-cξ[1])/rx)^2 + ((y-cξ[2])/ry)^2 -1
ξinit(xξ) = 0.5*(1-tanh(xξ*sqrt(W/(2*κ₀))))
#initvalξ = map((x,y)->(x<=5μm) ? 1.0 : (1.0-tanh(x-5μm*sqrt(W/(2.0*κ₀)))),grid)

initvalξ = map((x,y)->(x<=1) ? 1.0 : ξinit(xξ(x,y)),grid)
#initvalξ = map((x,y)-> ξinit(xξ(x,y)),grid)


vis1=GridVisualizer(Plotter=Plotter;
resolution=(500,400),
layout=(1,1), legend = :lt
)
scalarplot!(vis1,grid,initvalξ)


inivalLi[ϕₙ,:] .= 0.0
inivalLi[ζₙ,:] .= 1.0
inivalLi[ξₙ,:] = initvalξ 


problemLi = ODEProblem(sysLi,inivalLi,(0.0,10.0))

#=
solve(problemLi, KenCarp47(linsolve = KrylovJL_GMRES()),
    save_everystep = false);
using Symbolics
du0 = copy(inivalLi)
jac_sparsity = Symbolics.jacobian_sparsity((du, u) -> sysLi,
    du0, inivalLi)
tp=reshape(problemLi.u0,sysLi)
=#
tsol = solve(problemLi, Rodas5();
             adaptive = true,
             reltol = 1.0e-3,          # Tighter relative tolerance
             abstol = 1.0e-3,          # Tighter absolute tolerance
             initializealg = NoInit(),# Automatically choose the best initialization algorithm
            # force_dtmin = true,       # Force the minimum time step size
           #  dtmin = 1.0e-6,          # Set a reasonable minimum time step size
           #maxiters = 1000,           # Set a limit for maximum iterations
             #progress = true,          # Enable progress monitoring
             #progress_steps = 1      # Set progress steps
            ) 

#=
OrderedDict{String, UnionAll} with 5 entries:
  "QNDF2"                     => QNDF2
  "FBDF"                      => FBDF
  "Rosenbrock23 (Rosenbrock)" => Rosenbrock23
  "Implicit Euler"            => ImplicitEuler
  "Implicit Midpoint"         => ImplicitMidpoint
  "Rodas"                     => Rodas5
=#

tsol=reshape(tsol,sysLi)

vis1=GridVisualizer(Plotter=Plotter;
resolution=(1200,300),
layout=(1,3), legend = :lt
)



ξval = tsol.u[end][ξₙ,:]
ζval = tsol.u[end][ζₙ,:]
ϕval = tsol.u[end][ϕₙ,:]
scalarplot!(vis1[1, 1], grid, ξval)
scalarplot!(vis1[1, 2], grid, ζval)
scalarplot!(vis1[1, 3], grid, ϕval)
reveal(vis1)




#odesol = solve(problemfiber,Rosenbrock23())
#odesol = solve(problemfiber,QNDF2())
#odesol = solve(problemfiber,FBDF())
#odesol = solve(problemfiber,ImplicitEuler())

#=
function make_grid(;maxvolume=2)
	builder=SimplexGridBuilder(Generator=Triangulate)
	p00=point!(builder, 0,0)
	p10=point!(builder, 200,0.0)
	p11=point!(builder, 200,50)
	p01=point!(builder, 0,50)
	
	facetregion!(builder,1)
	facet!(builder, p00,p10)
	facetregion!(builder,2)
	facet!(builder, p10,p11)
	facetregion!(builder,3)
	facet!(builder, p11,p01)
	facetregion!(builder,4)
	facet!(builder,p00,p01)
		
	simplexgrid(builder,maxvolume=maxvolume)
end

grid2 = make_grid()

gridplot(grid2;Plotter=Plotter,resolution = (600, 250), legend =:rt)
=#
#@info "UMFPACK:"
#umf_sol = solve(SmoldSystem; inival = inivalSmold, method_linear = UMFPACKFactorization())
#=
sol = solve(SmoldSystem; 
method_linear = KrylovJL_BICGSTAB(),
precon_linear = UMFPACKFactorization(),
keepcurrent_linear =false,
)
=#


