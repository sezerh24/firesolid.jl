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
Plotter = GLMakie

#= M. Looyeh, P. Bettess, A. Gibson, A one-dimensional finite element simulation for
the fire-performance of GRP panels for offshore structures, Int. J. Numer. Methods
Heat Fluid Flow 7 (1997) 609–625.
=#


@unitfactors J K kg mol W m kJ °C s kW R cm

#= 
Full geometry
H_Fastener = 1.2cm
Width = 2cm 
Ny = 32
Nx = 32
FTW = 1cm
FBW = 0.5cm
X = collect(0:Width/Nx:Width)
Y = collect(0:H_Fastener/Ny:H_Fastener)
grid = simplexgrid(X, Y)
rect!(grid, [0.0, H_Fastener/4], [Width, H_Fastener]; region = 2, bregions = [1,2,3,2])
rect!(grid, [(Width-FTW)/2, 0.0], [(Width+FTW)/2, H_Fastener/4]; region = 3, bregions = [1, 1, 4, 1])
rect!(grid, [(Width-FBW)/2, H_Fastener/4], [(Width+FBW)/2, H_Fastener]; region = 3, bregions = [0, 4, 3, 4])
gridnew = subgrid(grid, [2,3])
=#

#Symetric Geometry
H_Fastener = 1.2cm
Width = 2cm/2
Ny = 64
Nx = 64
FTW = 1cm/2
FBW = 0.5cm/2
X = collect(0:Width/Nx:Width)
Y = collect(0:H_Fastener/Ny:H_Fastener)
grid = simplexgrid(X, Y)
rect!(grid, [0.0, H_Fastener/4], [Width, H_Fastener]; region = 2, bregions = [1,2,3,2])
rect!(grid, [(Width-FTW), 0.0], [Width, H_Fastener/4]; region = 3, bregions = [1, 1, 4, 1])
rect!(grid, [(Width-FBW), H_Fastener/4], [Width, H_Fastener]; region = 3, bregions = [0, 2, 3, 4])
gridnew = subgrid(grid, [2,3])


vis=GridVisualizer(Plotter=Plotter;
resolution=(600, 250),
layout=(1,1), legend = :rb, aspect = 1)
gridplot!(vis, gridnew, clear=true, show=true)
reveal(vis)


subgrid1 = subgrid(gridnew, [2])
subgrid2 = subgrid(gridnew, [3])

final_t = 1.0s
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
#Nx = 60;
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


kᵢ  = 1/(Vf/kf + (1-Vf)/kᵣ)
Cₚᵢ = 1/(Vf/Cₚf + (1-Vf)/Cₚr)

# if the fastener is 304 Stainless Steel (Austenitic)
kₛ = 16.2W/m/K
ρₛ = 8000kg/m^3
Cₚₛ = 500J/kg/K

#mᵢ = similar(X)
#ρΔ = similar(X)
#TΔ = similar(X)
ρₙ = 1;
Tₙ = 2; 
mₙ = 3; 
SpecInd = [ρₙ,Tₙ,mₙ]
function reactionT!(f,u,node,data)  

    if node.region == 2
        f[ρₙ] = A*u[ρₙ]*abs(exp((-Eₐ/R/u[Tₙ])))
        h     = Cₚᵢ*(u[Tₙ]-Tₐ)
        hg    = Cₚg*(u[Tₙ]-Tₐ)
	    f[Tₙ]  = A*u[ρₙ]*abs(exp((-Eₐ/R/u[Tₙ])))*(Q*10+h-hg)	
        f[mₙ]  = -A*u[ρₙ]*abs(exp((-Eₐ/R/u[Tₙ])))
    elseif node.region == 3
        f[Tₙ] = 0.0 
    end
    
end

function storageT!(f,u,node, data)
   
    if node.region == 2
        ρT    = u[ρₙ]*(1-Vf) + ρf*Vf
        f[ρₙ] = u[ρₙ]
        f[Tₙ] = ρT*u[Tₙ]*Cₚᵢ
    elseif node.region == 3
        f[Tₙ] = ρₛ*u[Tₙ]*Cₚₛ
    end
   
end

function fluxT!(f,u,edge,data)

    # convective and diffusive flux 
    #=
    velT        = project(edge,Cₚg*(u[mₙ,1] + u[mₙ,2])/2)
    f[mₙ]       = project(edge,u[mₙ,1]) 
    Bplus =  kᵢ*fbernoulli(velT/kᵢ)
    Bminus =  kᵢ*fbernoulli(velT/kᵢ)
    f[Tₙ] = Bminus * u[Tₙ, 1] - Bplus * u[Tₙ, 2]
    =#
    #println("edge region = ", edge.region)
    
 
   

    if edge.region == 2
        #println("Edge coordinates for region 2: ", edge.coord)
        #println("Validating mₙ and u[mₙ, 1]: ", mₙ, u[mₙ, 1])

        #println("edge region = ", edge.region)
	    # convective and diffusive flux 
        #f[mₙ]       = project(edge,u[mₙ,1])
        #println("edge region = ", edge.region)
       # velT        = project(edge,Cₚg*(u[mₙ,1] + u[mₙ,2])/2)

        velT   = Cₚg*(u[mₙ,1] + u[mₙ,2])/2
        f[mₙ]   = u[mₙ,1]
        
        Bplus =  kᵢ*fbernoulli(velT/kᵢ)
        Bminus =  kᵢ*fbernoulli(velT/kᵢ)
        f[Tₙ] = Bminus * u[Tₙ, 1] - Bplus * u[Tₙ, 2]

    elseif edge.region == 3
      #  println("edge region = ", edge.region)

        f[Tₙ] = kₛ*(u[Tₙ, 1] - u[Tₙ, 2])

    end
   

end

Tf(t) = Tₐ + (1100-(Tₐ-273.18))*(1-exp(-exp(0.71*log(t/124.8)))) - 100*(1-exp(-exp(0.71*log(t/124.8))))

function bconditionsT!(f,u,bnode,data)
	boundary_dirichlet!(f,u,bnode,species=Tₙ,region=1,value=Tf(bnode.time))
    boundary_dirichlet!(f,u,bnode,species=mₙ,region=3,value=0)
    #boundary_neumann!(f,u,bnode,species=mₙ,region=4,value=0)
    #boundary_neumann!(f,u,bnode,species=ρₙ,region=4,value=0)
end

data = (table = tabledata.Value)

sysfiber  =   VoronoiFVM.System(gridnew;data = data,
    flux       = fluxT!,
    storage    = storageT!,
    reaction   = reactionT!,
	bcondition = bconditionsT!,
    unknown_storage = :sparse,
    assembly = :cellwise)
# activate the species in the fuel electrode: Region: 1

enable_species!(sysfiber, SpecInd[Tₙ], [2,3])	
enable_species!(sysfiber, SpecInd[ρₙ], [2])
enable_species!(sysfiber, SpecInd[mₙ], [2])

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
problemfiber = ODEProblem(sysfiber,inivalfiber,(0,final_t*1000))
#=
using Symbolics
du0 = copy(inivalfiber)
jac_sparsity = Symbolics.jacobian_sparsity((du, u) -> sysfiber,
    du0, inivalfiber)

#f = ODEFunction(sysfiber; jac_prototype = float.(jac_sparsity))
=#
#=
tstep = 1.0
control = VoronoiFVM.NewtonControl()
control.verbose = false
control.Δt_min = tstep
control.Δt = tstep
control.Δt_grow = 1.1
control.Δt_max = 1
control.Δu_opt = 2.0
control.damp_initial = 0.5
tsol = solve(sysfiber;
                     method_linear = UMFPACKFactorization(),
                     inival=inivalfiber,
                     times = [0.0, 1000],
                     control = control,)
=#
#=
tsol=solve(problemfiber,QNDF2(),  

                                   force_dtmin=true,
                                   adaptive=true,
                                   reltol=1.0e-3,
                                   abstol=1.0e-3,
                                   #initializealg=NoInit(),
                                  #dtmin=1e-6,
                                  # force_dtmin = true,
                                   #progress=true,
                                   #progress_steps=1,

                                   )
                                   
=#
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


odesol = solve(problemfiber,Rosenbrock23())
#odesol = solve(problemfiber,QNDF2())
#odesol = solve(problemfiber,FBDF())
#odesol = solve(problemfiber,ImplicitEuler())
#tsol=reshape(tsol,sysfiber)
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
smoldsol=reshape(odesol,sysfiber)
concentrations = Dict{String, Vector{Float64}}()

tsol = smoldsol
#density = tsol.u[end][1, :]
Temp    = tsol.u[800][2, :]
#mass    = tsol.u[end][3, :]
density = view(tsol.u[800][1, :], subgrid1)
mass = view(tsol.u[800][3, :], subgrid1)

#=
ρₙ = 1;
Tₙ = 2; 
mₙ = 3; 
=#

#=
for (i, species) in enumerate(SpecInd)
    concentrations[species] = tsol.u[end][i, :] #/ FEtsol.u[end][i, 1]
end
=#
	# Initialize the plot
vis1=GridVisualizer(Plotter=Plotter;
resolution=(1500,400),
layout=(1,3), legend = :lt
)

scalarplot!(vis1[1,1], gridnew, Temp, label="Temperature")
scalarplot!(vis1[1,2],subgrid1, density, label="desnity")
scalarplot!(vis1[1,3], subgrid1, mass, label="mass")  

Tcell1 = [matrix[Tₙ, 1] for matrix in tsol.u] .-273.18
Tcell2 = [matrix[Tₙ, Int(floor(Nx/10))] for matrix in tsol.u] .-273.18
Tmid = [matrix[Tₙ, Int(floor(Nx/2))] for matrix in tsol.u] .-273.18
Tend = [matrix[Tₙ, Nx] for matrix in tsol.u] .-273.18

ρcell1 = [matrix[ρₙ, 1] for matrix in tsol.u]
ρcell2 = [matrix[ρₙ, Int(floor(Nx/10))] for matrix in tsol.u]
ρmid = [matrix[ρₙ, Int(floor(Nx/2))] for matrix in tsol.u]
ρend = [matrix[ρₙ, Nx] for matrix in tsol.u]

mdot = [matrix[mₙ, Int(floor(Nx/2))] for matrix in tsol.u]

plot(tsol.t,Tcell1, label="hot surface")
plot!(tsol.t,Tcell2)
plot(tsol.t,mdot)
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