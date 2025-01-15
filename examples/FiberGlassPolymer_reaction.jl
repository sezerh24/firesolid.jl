using VoronoiFVM
using ExtendableGrids
using ExtendableSparse
using GridVisualize
using LessUnitful
using LinearAlgebra
using LinearSolve
using StaticArrays
using DataFrames
using GLMakie
using Colors
using PrettyTables
using OrdinaryDiffEq
using MAT
Plotter = GLMakie


#= M. Looyeh, P. Bettess, A. Gibson, A one-dimensional finite element simulation for
the fire-performance of GRP panels for offshore structures, Int. J. Numer. Methods
Heat Fluid Flow 7 (1997) 609–625.
=#

@unitfactors J K kg mol W m kJ °C s kW R cm
final_t = 1000.0s
# Define the data as a DataFrame
tabledata = DataFrame(
    Property = [
        "Pre-exponential factor, A (1/sec)",
        "Specific heat (glass-fibre), (cp)fr (J/kg K)",
        "Specific heat (polyester resin), (cp)r (J/kg K)",
        "Specific heat (gas), cpg (J/kg K)",
        "Activation energy, EA (kJ/kg-mole)",
        "Thermal conductivity (glass-fibre, bulk), k fr (W/m K)",
        "Thermal conductivity (polyester resin), kr (W/m K)",
        "Heat of decomposition, Q (J/kg)",
        "Gas constant, R (kJ/kg-mole K)",
        "Ambient temperature, T° (°C)",
        "Volume fraction, Vf",
        "Density (initial), ro (kg/m3)",
        "Density (glass fibre), rfr (kg/m3)",
        "Density (polyester resin), rr (kg/m3)"
    ],
    Value = [1.0e3, 760.0, 1600.0, 2386.5, 0.5e5, 1.04, 0.20, 2.3446e5, 8.314, 20.0, 0.45, 1812.0, 2560.0, 1200.0]
)

# Print the table using PrettyTables
pretty_table(tabledata)  # Specify backend=:html for Jupyter notebooks
Nx = 60;
A = tabledata.Value[1]/s
Cₚf = tabledata.Value[2]J/kg/K 
Cₚr = 2*tabledata.Value[3]J/kg/K 
Cₚg = tabledata.Value[4]J/kg/K
Eₐ = tabledata.Value[5]J/kg/mol 
kf = tabledata.Value[6]W/m/K 
kᵣ = tabledata.Value[7]W/m/K
Q  = tabledata.Value[8]J/kg 
Tₐ = tabledata.Value[10]+273.18K 
Vf = tabledata.Value[11]
ρᵢ = tabledata.Value[12]kg/m^3
ρf = tabledata.Value[13]kg/m^3
ρr   = tabledata.Value[14]kg/m^3
L = 1.1cm 
h = L/(Nx-1); #grid size
X = 0:h:L
grid = simplexgrid(X)
gridplot(grid;Plotter=GLMakie,resolution = (600, 250), legend =:rt)

kᵢ  = 1/(Vf/kf + (1-Vf)/kᵣ)
Cₚᵢ = 1/(Vf/Cₚf + (1-Vf)/Cₚr)


mᵢ = similar(X)
ρΔ = similar(X)
TΔ = similar(X)

ρₙ = 1;
Tₙ = 2; 
mₙ = 3; 
SpecInd = [ρₙ,Tₙ,mₙ]
function reactionT!(f,u,node,data)  
    f[ρₙ] = A*u[ρₙ]*abs(exp((-Eₐ/R/u[Tₙ])))
    h     = Cₚᵢ*(u[Tₙ]-Tₐ)
    hg    = Cₚg*(u[Tₙ]-Tₐ)
	f[Tₙ]  = A*u[ρₙ]*abs(exp((-Eₐ/R/u[Tₙ])))*(Q*10+h-hg)	
    f[mₙ]  = -A*u[ρₙ]*abs(exp((-Eₐ/R/u[Tₙ])))
end

function storageT!(f,u,node, data)
     ρT   = u[ρₙ]*(1-Vf) + ρf*Vf
    f[ρₙ] = u[ρₙ]
    f[Tₙ] = ρT*u[Tₙ]*Cₚᵢ
   
   #=
    h     = Cₚᵢ*(u[Tₙ]-Tₐ)
    hg    = Cₚg*(u[Tₙ]-Tₐ)
    f[Tₙ] = ρT*u[Tₙ]*Cₚᵢ - f[ρₙ]*(Q*10+h-hg)
   =#
end

function fluxT!(f,u,edge,data)

	# convective and diffusive flux 
    velT   = project(edge,Cₚg*(u[mₙ,1] + u[mₙ,2])/2)
    f[mₙ]  = project(edge,u[mₙ,1]) 
    #velT   = Cₚg*(u[mₙ,1] + u[mₙ,2])/2
    #f[mₙ]   = u[mₙ,1]
    Bplus  =  kᵢ*fbernoulli(velT/kᵢ)
    Bminus =  kᵢ*fbernoulli(velT/kᵢ)
    f[Tₙ]   = Bminus * u[Tₙ, 1] - Bplus * u[Tₙ, 2]

end
#=
function fluxT!(f,u,edge,data)

	# Flux for ρg
    #=
    nk = edge.node[1]
    nl = edge.node[2]
    Δx = edge.coord[nl]-edge.coord[nk]
    #println("Δx  = ", Δx)
    ρΔ[nk] = (u[ρₙ,1] + u[ρₙ,2])/2
    TΔ[nk] = (u[Tₙ,1] + u[Tₙ,2])/2
  
     mm = 0;
    for i = nk:Nx-1
    #for i = 1:nk
        #println("i = ", i)
        mm +=  A*ρΔ[i]*abs(exp((-Eₐ/R/TΔ[i])))*Δx
    end
    =#
   # mm =  A*ρΔ[nk]*abs(exp((-Eₐ/R/TΔ[nk])))*Δx

   # mm =  A*u[rho]*abs(exp((-Eₐ/R/TΔ[nk])))*Δx

    
    #= 
    m = 0.0
    for i = 1:nk
        m += -Cₚg*(A*u[ρₙ,1]*exp(-Eₐ/R/u[Tₙ,1])*(X[nk]-L) + A*u[ρₙ,2]*exp(-Eₐ/R/u[Tₙ,2])*(X[nl]-L))/2
    end
    =#
	# convective and diffusive flux 
    velT        = project(edge,Cₚg*(u[mₙ,1] + u[mₙ,2])/2)
    #velT = -mm*Cₚg
    #=convection  = velT>0 ? velT*u[Tₙ,1] : velT*u[Tₙ,2] 
	f[Tₙ] 		= kᵢ*(u[Tₙ,1]-u[Tₙ,2]) + convection =#
    f[mₙ]       = project(edge,u[mₙ,1]) 
    
    #=
    t = edge.time
    println("edge time = ", t)
    =#
    Bplus =  kᵢ*fbernoulli(velT/kᵢ)
    Bminus =  kᵢ*fbernoulli(velT/kᵢ)
    f[Tₙ] = Bminus * u[Tₙ, 1] - Bplus * u[Tₙ, 2]

end
=#

Tf(t) = Tₐ + (1100-(Tₐ-273.18))*(1-exp(-exp(0.71*log(t/124.8)))) - 100*(1-exp(-exp(0.71*log(t/124.8))))

function bconditionsT!(f,u,bnode,data)

	

    # Here we access the boundary node and give the time dependent bcondition
   # v=ramp(bnode.time,dt=(0,final_t),du=(Tf(bnode.time), Tf(bnode.time)))

  # v=ramp(bnode.time,dt=(0,final_t),du=(Tf(bnode.time), Tf(bnode.time)))
   # v=ramp(bnode.time,dt=(0,final_t),du=(1000.0, Tf(bnode.time)))
   
  # println("boundary v = ", v)
   #=
    Rb = 1-exp(-exp(0.71*log(t/124.8)))
    Tf = Tₐ + (1100-Tₐ)*Rb - 100*Rb =#
#    v=ramp(bnode.time,dt=(0,tg),du=(value, 0.0))
	boundary_dirichlet!(f,u,bnode,species=Tₙ,region=1,value=Tf(bnode.time))
    boundary_dirichlet!(f,u,bnode,species=mₙ,region=2,value=0)
end

data = (table = tabledata.Value)

sysfiber  =   VoronoiFVM.System(grid;data = data,
    flux       = fluxT!,
    storage    = storageT!,
    reaction   = reactionT!,
	bcondition = bconditionsT!,
    unknown_storage = :sparse,
    assembly = :edgewise)
# activate the species in the fuel electrode: Region: 1

for i in SpecInd
    enable_species!(sysfiber, SpecInd[i], [1])	
end

inivalfiber=unknowns(sysfiber)
inivalfiber[Tₙ,:] .= Tₐ
inivalfiber[ρₙ,:] .= ρr
inivalfiber[mₙ,:] .= 0.0

#@info "UMFPACK:"
#umf_sol = solve(SmoldSystem; inival = inivalSmold, method_linear = UMFPACKFactorization())
#=
sol = solve(sysfiber; 
method_linear = KrylovJL_BICGSTAB(),
precon_linear = UMFPACKFactorization(),
keepcurrent_linear =false,
)=#

#dt = 1
problemfiber = ODEProblem(sysfiber,inivalfiber,(0,final_t))
#=
using Symbolics
du0 = copy(inivalfiber)
jac_sparsity = Symbolics.jacobian_sparsity((du, u) -> sysfiber,
    du0, inivalfiber)

#f = ODEFunction(sysfiber; jac_prototype = float.(jac_sparsity))
=#
tsol=solve(problemfiber,ImplicitEuler(),  

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
                                   

#=
 NOTE: Rodas5P works very well even if a dtmin is not given. However, implicitEuler requires a dtmin depending on the problem being solved... 
Rodas5Pr is working, however, it does not adapt the time as good ad Rodas5P and also Rodas5Pr is not working when a dtmin is not defined. an
appropriate dtmin should be defined otherwise the method does not work. 
ImplicitEuler works with dtmin = 10s. however, the other solvers does not work with such a high time step... 



OrderedDict{String, UnionAll} with 5 entries:
  "QNDF2"                     => QNDF2
  "FBDF"                      => FBDF
  "Rosenbrock23 (Rosenbrock)" => Rosenbrock23
  "Implicit Euler"            => ImplicitEuler
  "Implicit Midpoint"         => ImplicitMidpoint
=#


#odesol = solve(problemfiber,Rosenbrock23())
#odesol = solve(problemfiber,QNDF2())
#odesol = solve(problemfiber,FBDF())
#odesol = solve(problemfiber,ImplicitEuler())
tsol=reshape(tsol,sysfiber)
#=
Sol = solve(
    smold_ode_problem,
		ImplicitEuler(),
		adaptive=true,
        #dt = 5,
        reltol=1.0e-3,
		abstol=1.0e-3,
		initializealg=NoInit()
	)	
smoldsol=reshape(Sol,SmoldSystem)
=#
	# Initialize the plot
vis1=GridVisualizer(Plotter=Plotter;
resolution=(400,300),
layout=(1,3), legend = :lt
)
#scalarplot!(vis1[1,1], grid, tsol.u[end][2,:], label="Temperature")
Tcell1 = [matrix[Tₙ, 1] for matrix in tsol.u] .-273.18
Tcell2 = [matrix[Tₙ, Int(floor(Nx/10))] for matrix in tsol.u] .-273.18
Tmid = [matrix[Tₙ, Int(floor(Nx/2))] for matrix in tsol.u] .-273.18
Tend = [matrix[Tₙ, Nx] for matrix in tsol.u] .-273.18

ρcell1 = [matrix[ρₙ, 1] for matrix in tsol.u]
ρcell2 = [matrix[ρₙ, Int(floor(Nx/10))] for matrix in tsol.u]
ρmid = [matrix[ρₙ, Int(floor(Nx/2))] for matrix in tsol.u]
ρend = [matrix[ρₙ, Nx] for matrix in tsol.u]

mdot = [matrix[mₙ, Int(floor(Nx/2))] for matrix in tsol.u]

scalarplot!(vis1[1,1], tsol.t,Tcell1, label="hot surface", markershape=:cross, color=:black, markevery=1, linestyle=:dash, clear = false,
xlabel="time / s", ylabel="Temperature")
scalarplot!(vis1[1,1], tsol.t,Tcell2 , label="x/L = 0.1", markershape=:cross, color=:blue, markevery=1, linestyle=:dash, clear = false,
xlabel="time / s", ylabel="Temperature")
scalarplot!(vis1[1,1], tsol.t,Tmid , label="x/L = 0.5", markershape=:cross, color=:red, markevery=1, linestyle=:dash, clear = false,
xlabel="time / s", ylabel="Temperature")
scalarplot!(vis1[1,1], tsol.t,Tend , label="cold surface", markershape=:cross, color=:green, markevery=1, linestyle=:dash, clear = false,
xlabel="time / s", ylabel="Temperature")

# Function to generate a list of colors
function generate_colors(n::Int)
    colors = [:red, :green, :blue, :purple, :cyan, :magenta, :brown, :black]
    return colors
end

# Load data from .mat file
ExpData = matread("ExpData.mat")

# Extract data arrays
x_sim_HotSide = ExpData["x_sim_HotSide"]
y_sim_HotSide = ExpData["y_sim_HotSide"]
x_Sim_xL1_10 = ExpData["x_Sim_xL1_10"]
y_Sim_xL1_10 = ExpData["y_Sim_xL1_10"]
x_ex_xL1_10 = ExpData["x_ex_xL1_10"]
y_ex_xL1_10 = ExpData["y_ex_xL1_10"]
x_Sim_xL5_10 = ExpData["x_Sim_xL5_10"]
y_Sim_xL5_10 = ExpData["y_Sim_xL5_10"]
x_ex_xL5_10 = ExpData["x_ex_xL5_10"]
y_ex_xL5_10 = ExpData["y_ex_xL5_10"]
x_sim_ColdSide = ExpData["x_sim_ColdSide"]
y_sim_ColdSide = ExpData["y_sim_ColdSide"]
x_ex_ColdSide = ExpData["x_ex_ColdSide"]
y_ex_ColdSide = ExpData["y_ex_ColdSide"]

scalarplot!(vis1[1,1], x_sim_ColdSide,y_sim_ColdSide,label="sim-cold surface", markershape=:cross, color=:purple, markevery=5, linestyle=:dash, clear = false,
xlabel="time / s", ylabel="Temperature")
scalarplot!(vis1[1,1], x_ex_ColdSide,y_ex_ColdSide,label = "ex-cold surface", markershape=:cross, color=:magenta, markevery=5, linestyle=:dash, clear = false,
xlabel="time / s", ylabel="Temperature")

scalarplot!(vis1[1,1], x_Sim_xL5_10,y_Sim_xL5_10,label="sim-x/L=0.5", markershape=:cross, color=:purple, markevery=5, linestyle=:dash, clear = false,
xlabel="time / s", ylabel="Temperature")
scalarplot!(vis1[1,1], x_ex_xL5_10,y_ex_xL5_10,label = "ex-x/L=0.5", markershape=:cross, color=:magenta, markevery=5, linestyle=:dash, clear = false,
xlabel="time / s", ylabel="Temperature")

scalarplot!(vis1[1,1], x_Sim_xL1_10,y_Sim_xL1_10,label="sim-x/L=0.1", markershape=:cross, color=:purple, markevery=5, linestyle=:dash, clear = false,
xlabel="time / s", ylabel="Temperature")
scalarplot!(vis1[1,1], x_ex_xL1_10,y_ex_xL1_10,label = "ex-x/L=0.1", markershape=:cross, color=:magenta, markevery=5, linestyle=:dash, clear = false,
xlabel="time / s", ylabel="Temperature")

scalarplot!(vis1[1,2], tsol.t,ρcell1 , label="hot surface", markershape=:cross, color=:black, markevery=1, linestyle=:dash, clear = false,
xlabel="time / s", ylabel="ρ")
scalarplot!(vis1[1,2], tsol.t,ρcell2 , label="x/L = 0.1", markershape=:cross, color=:blue, markevery=1, linestyle=:dash, clear = false,
xlabel="time / s", ylabel="ρ")
scalarplot!(vis1[1,2], tsol.t,ρmid , label="x/L = 0.5", markershape=:cross, color=:red, markevery=1, linestyle=:dash, clear = false,
xlabel="time / s", ylabel="ρ")
scalarplot!(vis1[1,2], tsol.t,ρend , label="cold surface", markershape=:cross, color=:green, markevery=1, linestyle=:dash, clear = false,
xlabel="time / s", ylabel="ρ")


scalarplot!(vis1[1,3], tsol.t,mdot , label="ntot", markershape=:cross, color=:black, markevery=1, linestyle=:dash, clear = false,
xlabel="time / s", ylabel="mdot")

reveal(vis1)
#func = similar(Tvec)
#func .= A.*abs.(ρvec.*exp.((-Eₐ/A/R./Tvec)))