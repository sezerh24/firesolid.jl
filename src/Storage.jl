
function Fstorage_Array!(f, u, node,data)
    # storage term for the gaseous species inside the fuel electrode
    nn  = node.index # accessing the grid nodes
    # n is accessing number of species... 
    inds = data.E["nSp"]
    ϵ   =data.E["ϵ"][nn]
    f[1:inds] = u[1:inds]*ϵ
end  

function Fstorage!(f, u, node,data)
    # storage term for the gaseous species inside the fuel electrode
    #nn  = node.index # accessing the grid nodes
    # n is accessing number of species... 
    inds = data.E["nSp"]
    ϵ   =data.E["ϵ"]
    f[1:inds] = u[1:inds]*ϵ
end  

function Fstorage_Epot_Array!(f, u, node,data)

    SI = data.SpecIndE;
    E = data.E
    
    # storage term for the gaseous species inside the fuel electrode
    nn  = node.index # accessing the grid nodes
    inds = E["nSp"]

    ϵ   = E["ϵ"][nn];
    Cdl = E["Cdl"][nn];
    f[1:inds] = u[1:inds]*ϵ

    # here we add the double capacitor for electrical potential in anode 
    f[SI["ϕaₑ"]] = Cdl*(u[SI["ϕaₑ"]] - u[SI["ϕᵢ"]])
    # here we add the double capacitor for ionic potential in anode 
    f[SI["ϕᵢ"]] =  Cdl*(u[SI["ϕᵢ"]] - u[SI["ϕaₑ"]])

end


function Fstorage_Epot!(f, u, node,data)

    SI = data.SpecIndE;
    E  = data.E
    
    # storage term for the gaseous species inside the fuel electrode
    #nn  = node.index # accessing the grid nodes
    inds = E["nSp"]
    ϵ    = E["ϵ"];
    Cdl  = E["Cdl"];
    f[1:inds] = u[1:inds]*ϵ
    # here we add the double capacitor for electrical potential in anode 
    f[SI["ϕₑ"]] = Cdl*(u[SI["ϕₑ"]] - u[SI["ϕᵢ"]])
    # here we add the double capacitor for ionic potential in anode 
    f[SI["ϕᵢ"]] =  Cdl*(u[SI["ϕᵢ"]] - u[SI["ϕₑ"]])

end


