
function Freaction_Epot_Array!(f,u,node,data)

    SI = data.SpecIndE;
    E  = data.E; 
    R  = data.R
    T  = data.T;
    F  = data.F;

    nn      = node.index
    lTPB    = E["lTPB"][nn]
    H2inf   = E["ImF"][SI["H2"]]   # should be mole/m^3
    H2Oinf  = E["ImF"][SI["H2O"]]  # should be mole/m^3

    iₒ = 31.4*(R*T*u[SI["H2"]])^(-0.03)*(R*T*u[SI["H2O"]])^(0.4)*exp(-18300/T)
    ηₐ = u[SI["ϕₑ"]] - u[SI["ϕᵢ"]]
         - (R*T/(2*F))*
         log((H2inf/u[SI["H2"]])*(u[SI["H2O"]]/H2Oinf))
		 
    iFBV = iₒ*lTPB*(exp((F/(R*T))*ηₐ)-exp((-F/(R*T))*ηₐ)) 
    f[SI["ϕᵢ"]] = -iFBV
	f[SI["ϕaₑ"]] = iFBV
	f[SI["H2"]] = iFBV/(2*F)  # the source term for the hydrogen 
	f[SI["H2O"]] = -iFBV/(2*F) # the source term for the steam
end


function Freaction_Epot!(f,u,node,data)

    SI = data.SpecIndE;
    E  = data.E; 
    R  = data.R
    T  = data.T;
    F  = data.F;

    lTPB    = E["lTPB"]
    H2inf   = E["ImF"][SI["H2"]]   # should be mole/m^3
    H2Oinf  = E["ImF"][SI["H2O"]]  # should be mole/m^3

    iₒ = 31.4*(R*T*u[SI["H2"]]+0im)^(-0.03)*(R*T*u[SI["H2O"]]+0im)^(0.4)*exp(-18300/T)
    ηₐ = u[SI["ϕₑ"]] - u[SI["ϕᵢ"]]
         - (R*T/(2*F))*
         log(Complex((H2inf/u[SI["H2"]])*(u[SI["H2O"]]/H2Oinf)))
		 
    iFBV = abs(iₒ*lTPB*(exp((F/(R*T))*ηₐ)-exp((-F/(R*T))*ηₐ)))
    f[SI["ϕᵢ"]] = -iFBV
	f[SI["ϕₑ"]] = iFBV
	f[SI["H2"]] = iFBV/(2*F)  # the source term for the hydrogen 
	f[SI["H2O"]] = -iFBV/(2*F) # the source term for the steam
end



function Fsource!(f,node,data)
    SI = data.SpecIndE;
    iF = data.iF
    nn  = node.index
    f[SI["H2"]]  = -iF[nn]  # the source term for the hydrogen 
    f[SI["H2O"]] = iF[nn] # the source term for the steam  
end

