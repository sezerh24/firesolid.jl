# Step 1: Define the Reaction data type
mutable struct Reaction
    species_id::Int
    reaction_func::Function
    args_ids::Vector{Int}
end

mutable struct Storage
    species_id::Int
    storage_func::Function
    storage_const
    args_ids::Vector{Int}
end


mutable struct StorageC
    species_id::Int
    storage_const #::Vector{Float64}  # Array of Float64 with 2 dimensions
end

mutable struct Source
    species_id::Int
    source_const #::Vector{Float64}  # Array of Float64 with 2 dimensions
end

mutable struct Fluxses
    species_id::Int
    diffusion_func::Function
    args_ids::Vector{Int}
end

mutable struct DiffusionC
    species_id::Int
    diffusion_const #::Vector{Float64}  # Array of Float64 with 2 dimensions
end

mutable struct ConvectionF
    species_id::Int
    convection_func::Function
    args_ids::Vector{Int}
end

mutable struct ConvectionC
    species_id::Int
    convection_const #::Vector{Float64}  # Array of Float64 with 2 dimensions
end



# Step 4: Define the reaction function that handles any number of arguments
function reaction!(f, u, node, data)

    reactions = data.RCT

    for i in eachindex(reactions)
        k = reactions[i].species_id
        kl = reactions[i].args_ids
        # Collect the arguments dynamically
        args = [u[arg_id] for arg_id in kl]
        # Call the reaction function with the collected arguments
        f[k] = reactions[i].reaction_func(args...)
    end
end


# Step 4: Define the reaction function that handles any number of arguments
#=
function storage!(f, u, node, data)

    nk =node.index[1]
    

    StF = data.StF
    StC = data.StC

    for i in eachindex(StF)
        k = StF[i].species_id
        kl = StF[i].args_ids
        # Collect the arguments dynamically
        args = [u[arg_id] for arg_id in kl]
        # Call the reaction function with the collected arguments
        f[k] = u[k]*StF[i].storage_func(args...)
    end

    for i in eachindex(StC)
        k = StC[i].species_id
        f[k] = u[i]*StC[i].storage_const[nk]
    end
end
=#

function storage!(f, u, node, data)
    StF = data.StF

    for i in eachindex(StF)
        k = StF[i].species_id
        kf = StF[i].args_ids
        # Collect the arguments dynamically
        args = [u[arg_id] for arg_id in kf]
        # Call the reaction function with the collected arguments
        f[k] = u[k]*(StF[i].storage_func(args...)+StF[i].storage_const)
    end

    #=
    for i in eachindex(StC)
        k = StC[i].species_id
        f[k] = u[i]*StC[i].storage_const
    end
    =#
end


function source!(f, node, data)
    Scr = data.Scr
    for i in eachindex(Scr)
        k = Scr[i].species_id
        f[k] = Scr[i].storage_const
    end
end

function source!(f, node, data)
    Scr = data.Scr
    nk = node.index[1]
    for i in eachindex(Scr)
        k = Scr[i].species_id
        f[k] = Scr[i].storage_const[nk]
    end
end



function Fluxfunc!(f, u, edge, data)

    nk = edge.index[1]
    nl = edge.index[2]

    DF = data.DF
    DC = data.DC
    CF = data.CF  # velocity is a function of the other variables 
    CC = data.CF  # velocity is a constant value
    CD = data.DF  # velocity is calculated from Darcy's flow 
    gasspec = data.gasspec # gas species 


    Pgrad = 0.0; 
    for i in eachindex(gasspec)
        Pgrad = Pgrad + (R*T)*(u[i, 1] - u[i, 2])  
    end

    for i in eachindex(DFC)
        k = DFC[i].species_id
        kd = DFC[i].Dargs_ids
        kc = DFC[i].Cargs_ids
        # Collect the arguments dynamically
        Dargs1 = [u[arg_id,1] for arg_id in kd]
        Dargs2 = [u[arg_id,2] for arg_id in kd]

        Cargs1 = [u[arg_id,1] for arg_id in kc]
        Cargs2 = [u[arg_id,2] for arg_id in kc]

        # Arithmetic mean of the diffusion
        D = (DFC[i].diffusion_func(Dargs1...) + DFC[i].diffusion_func(Dargs2...))/2
        V = (DFC[i].convection_func(Cargs1...) + DFC[i].convection_func(Cargs2...))/2

        Bplus = fbernoulli(V)
        Bminus = fbernoulli(V)
        # Call the reaction function with the collected arguments
        f[k] = D*(Bminus*u[k,1] - Bplus*u[k,2]) 
    end


    for i in eachindex(DC)
        k = DC[i].species_id
        D = (DC[i].diffusion_const[nk] + DC[i].diffusion_const[nl])/2
        f[k] = D*(u[k,1]-u[k,2])
    end

end
 






function Flux_Ficks_Darcy!(f, u, edge, data)

    nk1 = edge.node[1]
    nl1 = edge.node[2]

    E  = data.E; 
    R  = data.R
    T  = data.T;
    P = data.P;

     # Here calculate the multiComponent Diffusion 
     ϵ   = (E["ϵ"][nk1] + E["ϵ"][nl1])/2
     τ   = (E["τ"][nk1] + E["τ"][nl1])/2
     ΦPo = (E["pore_dia"][nk1] + E["pore_dia"][nl1])/2
     ΦPa = (E["part_dia"][nk1] + E["part_dia"][nl1])/2
     pm = Properties(ϵ,τ,ΦPo,ΦPa)
     # D_ij_e = (pm.ϵ/pm.τ) .* D_ij(E["tData"], T, P, E["tObj"].molwt)
     D_ije = D_ij(E["tData"], T, P, E["tObj"].molwt)
    for i = 1:E["nSp"]
        E["mF"][i] = (u[i,1]+u[i,2])/2.0;
    end

    totalMF = sum(E["mF"]);

    for i = 1:E["nSp"]
        E["mF"][i] = E["mF"][i]/totalMF;
    end

     # Knudsen Diffusion 
    # DK = (pm.ϵ/pm.τ) * D_Kn(E["tObj"].molwt, pm, T)
     DK =  D_Kn(E["tObj"].molwt, pm, T)
    # Here we calculate the molecular diffusion 
    # Initialization
    suma = 0.0
    Dm = zeros(E["nSp"])
    Mm = 0.0
    alpha = zeros(E["nSp"])
    Deff = zeros(E["nSp"])

    # First loop: Calculate Molecular diffusion of each specie
    for i in 1:E["nSp"]
        suma = 0.0
        for j in 1:E["nSp"]
            if i != j
                suma += E["mF"][j]/ (D_ije[i,j])
            end
        end
        Dm[i] = (1.0 - E["mF"][i]) / (suma)
    end

    # Second loop: Calculate average molecular weight in the mixture 
    for i in 1:E["nSp"]
        Mm += E["mF"][i] * E["tObj"].molwt[i]
    end

    # Third loop: Calculate the effective diffusion of the specie ith in the mixture 
    for i in 1:E["nSp"]
        alpha[i] = 1.0 - sqrt(E["tObj"].molwt[i] / (Mm))
        Deff[i] = (pm.ϵ/pm.τ) * (1.0 / ((1.0 - alpha[i] * E["mF"][i]) / (Dm[i]) + 1.0 / (DK[i])))
    end


    μ = viscosity(E["tData"],T,E["tObj"].molwt,E["mF"])
    # Calculate ∇Pₜ
    Pₜ = 0.0; 
    for i = 1:E["nSp"]
        Pₜ = Pₜ + (R*T)*(u[i, 1] - u[i, 2])  
    end

    for i = 1:E["nSp"]
                
            #vh = (K/(μ*DK[i]))*Pₜ
            vh = (K/(μ))*Pₜ
            Bplus = fbernoulli(vh)
            Bminus = fbernoulli(-vh)
            f[i]= Deff[i]*((Bminus*u[i,1])-Bplus*u[i,2])
 
    end  

end