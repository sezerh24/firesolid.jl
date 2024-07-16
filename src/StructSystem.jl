


mutable struct StorageStruct
    SpeciesId::Int
    Func::Function
    Const
    ArgsIds::Vector{Int}
end

mutable struct ReactionStruct
    SpeciesId::Int
    Func::Function
    ArgsIds::Vector{Int}
end

mutable struct SourceStruct
    SpeciesId::Int
    SourceConst #::Vector{Float64}  # Array of Float64 with 2 dimensions
end

mutable struct FluxStruct
    SpeciesId::Int
    FuncDiff::Function
    ArgsIdsDiff::Vector{Int}
    ConstDiff
    FuncV::Function
    ArgsIdsV
    ConstV
    DarcyId
end



function storage!(f, u, node, data)
    st = data.st
    nk = node.index[1]
    for i in eachindex(st)
        k = st[i].SpeciesId
        kf = st[i].ArgsIds
        # Collect the arguments dynamically
        args = [u[arg_id] for arg_id in kf]
        # Call the reaction function with the collected arguments
        Const = length(st[i].Const)>1 ? st[i].Const[nk] : st[i].Const
        f[k] = u[k]*(st[i].Func(args...)+Const)
    end
end

function reaction!(f, u, node, data)
    reactions = data.rxn
    for i in eachindex(reactions)
        k = reactions[i].SpeciesId
        kl = reactions[i].ArgsIds
        # Collect the arguments dynamically
        args = [u[arg_id] for arg_id in kl]
        # Call the reaction function with the collected arguments
        f[k] = reactions[i].Func(args...)
    end
end

function source!(f, node, data)
    scr = data.scr
    nk = node.index[1]
    for i in eachindex(scr)
        k = scr[i].SpeciesId
        Const = length(scr[i].Const)>1 ? scr[i].Const[nk] :
         scr[i].Const
        f[k] = Const
    end
end



function Fluxfunc!(f, u, edge, data)
    
    #=
    nk  = edge.index[1]
    nl  = edge.index[2]
    =#

    fst = data.FluxStruct
    gs  = data.gasspec # gas species 
    

    if MDiff == 1

        tData = data.tData
        tObj  = data.tObj
        pm    = data.pm
        P     = data.P
        T     = data.T
        R     = data.R

        for (i,j) in pairs(gs)
            SpmF[i] = (u[j,1]+u[j,2])/2.0;
        end
        totalMF = sum(SpmF);
        SpmF   .= SpmF./totalMF

        Deff, μ = MultiComponentFickeanDiff(pm, P,T,SpmF,tData,tObj)
        

        for (i,j) in pairs(gs)
            fst[j].ConstDiff = Deff[i]
        end
    end

    if Darcy == 1

        T     = data.T
        R     = data.R
        Permiability = data.Permiability

        Pgrad = 0.0; 
        for i in gs
            Pgrad = Pgrad + (R*T)*(u[i, 1] - u[i, 2])  
        end

        DarcyVel = -(Permiability/μ)*Pgrad 

        for i in (gs)
            fst[i].ConstV = DarcyVel
        end

    end

    for i in eachindex(fst)

        k  = fst[i].SpeciesId
        kd = fst[i].ArgsIdsDiff
        kc = fst[i].ArgsIdsV
        # Collect the arguments dynamically
        Dargs1 = [u[arg_id,1] for arg_id in kd]
        Dargs2 = [u[arg_id,2] for arg_id in kd]

        Cargs1 = [u[arg_id,1] for arg_id in kc]
        Cargs2 = [u[arg_id,2] for arg_id in kc]

        # Arithmetic mean of the diffusion
        #=
        Dconst = length(fst[i].ConstDiff)>1 ? (fst[i].ConstDiff[nk] + fst[i].ConstDiff[nl])/2 :
         fst[i].ConstDiff
         Vconst = length(fst[i].ConstV)>1 ? (fst[i].ConstV[nk] + fst[i].ConstV[nl])/2 :
         fst[i].ConstV
         =#
        D = fst[i].ConstDiff+(fst[i].FuncDiff(Dargs1...) + fst[i].FuncDiff(Dargs2...))/2
        V = fst[i].ConstV + (fst[i].FuncV(Cargs1...) + fst[i].FuncV(Cargs2...))/2
        Bplus = fbernoulli(V)
        Bminus = fbernoulli(V)
        # Call the reaction function with the collected arguments
        f[k] = D*(Bminus*u[k,1] - Bplus*u[k,2]) 
    end

end




function MultiComponentFickeanDiff(pm, P,T,SpmF,tData,tObj)
     # Binary diffusion/Molecular Diffusion
     D_ije = D_ij(tData, T, P, tObj.molwt)
    # Knudsen Diffusion 
     DK =  D_Kn(tObj.molwt, pm, T)
    # Here we calculate the molecular diffusion 
    # Initialization
    suma = 0.0
    Dm = zeros(SpmF)
    Mm = 0.0
    alpha = zeros(SpmF)
    Deff = zeros(SpmF)
    # First loop: Calculate Molecular diffusion of each specie
    for i in eachindex(SpmF)
        suma = 0.0
        for j in eachindex(SpmF)
            if i != j
                suma += SpmF[j]/ (D_ije[i,j])
            end
        end
        Dm[i] = (1.0 - SpmF[i]) / (suma)
    end

    # Second loop: Calculate average molecular weight in the mixture 
    for i in eachindex(SpmF)
        Mm += SpmF[i] * tObj.molwt[i]
    end

    # Third loop: Calculate the effective diffusion of the specie ith in the mixture 
    for i in eachindex(SpmF)
        alpha[i] = 1.0 - sqrt(tObj.molwt[i] / (Mm))
        Deff[i] = (pm.ϵ/pm.τ) * (1.0 / ((1.0 - alpha[i] * SpmF[i]) / (Dm[i]) + 1.0 / (DK[i])))
    end

    μ = viscosity(tData,T,tObj.molwt,SpmF)
    # Calculate ∇Pₜ
    return Deff, μ
end


function DGM_Diffusion(DK::Array{Float64,1}, Dij::Matrix{Float64}, Yᵢ)
    n = size(DK, 1)
    H = Matrix{Float64}(undef, n, n)
    for k=1:n
        sum = 0.0
        for l=1:n
            H[k,l] = -Yᵢ[k]/Dij[k,l]
            if k != l
                sum += Yᵢ[l]/Dij[k,l]
            end
        end
        H[k,k] = (1/DK[k]) + sum       
    end    
    H .= inv(H)
    return H
end

function MultiComponentDGM(pm, P, T,SpmF, tData, tObj)

    # Binary diffusion/Molecular Diffusion
    D_ij_e = (pm.ϵ/pm.τ)*D_ij(tData, T, P, tObj.molwt)
    # Knudsen Diffusion 
    DK =  (pm.ϵ/pm.τ) * D_Kn(tObj.molwt, pm, T)

    # Knudsen Diffusion 
    DK = D_Kn(E["tObj"].molwt, pm, T)
    # Dusty Gas Diffusion 
    DGM = DGM_Diffusion(DK, D_ij_e, SpmF)
    # Mixture viscosity
    μ = viscosity(E["tData"],T,E["tObj"].molwt,E["mF"]) 
   # end of diffusion calculations

    return DGM, μ
end



