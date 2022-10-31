function [ Psuc , Tsuc] = Sinlet(p_suc , T_suc , INport_Amax, INport_Amin, V_comp, MM_g , n_van,rpm,c,toll_d,c_v,aQ,bQ,cQ,pipe,cpitch,ct,lenght,D_i,D_f,roughness,mu_g)
% This function set up the iterative approach useful for the evaluation of
%the real temperature and pressures in the inlet of the air-end section 
%or the true beggining of the compression/expansion process
%
%
% INPUT
% P_suc: Suction pressure[Pa] at the begginning of the inlet duct
% T_suc: Inlet temperature [K] in correspondance of the inlet duct entrance
% INport_Amax: Inlet port maximum passage area  [m^2]
% INport_Amin: Inlet port minimum passage area  [m^2]
% V_comp(1): First closed cell volume
% MM_g: Molecular mass of the gas [kg/kmol]
% Rgas: Specific gas costant [J/kgK]
% n_van [-]           : number of vanes
% rpm [rpm]           : shaft angular speed
% c [-]               : geometry index
%   1                     : circular stator
%   2                     : elliptical stator
% % toll_d [-]          : discretization tolerance
% 
% OUTPUT
% P_suc: Real suction pressure[Pa] at the begginning of the compression
% process (inlet of the air-end)
% T_suc: Real inlet temperature [K] in correspondance of the air end inlet
% (beggining of the compression process)
% 
% 
% HISTORY:  Gianoncelli_Genoni: creation of the file, see thesis for
% further information

 %% DEFINITIONS %%
 port_type="inlet";
 R_g        = SX_Constant({'UniGasConstant'})/MM_g;          % specific gas constant [J/ kg K];
 gamma= (c_v + R_g)/c_v;                                     % 
 m_gas_guess = p_suc*V_comp(1)/(R_g*T_suc)*c*n_van*rpm/60;   % portata massica che esce dal modello di andre 
 alfa=0.00001;                                               % loop variable
 Loopinlet=1;
 toll_d=1e-6;
 fb=1;
 coeff=[aQ,bQ,cQ];
 mastseries=[0];
 counter=0;

 while Loopinlet

[p_f,T_f,rho_f,delta_p] = Concentrated_Losses(fb,coeff,m_gas_guess,p_suc,T_suc,MM_g,gamma);  %intake valve pressure drop

[p_f,T_f,rho_f,delta_p] = Distributed_Losses(fb,pipe,cpitch,ct,lenght,D_i,D_f,roughness,p_f,T_f,m_gas_guess,MM_g,gamma,mu_g);

[m_port_inf , Uinf , Pinf ,Tinf, Pout , Tout] = PortModel(p_f , T_f , INport_Amax , INport_Amin, gamma, R_g,  m_gas_guess, port_type);



mast = Pout*V_comp(1)/(R_g*Tout)*c*n_van*rpm/60;
mastseries=[mastseries , mast];

err= abs(mast-m_gas_guess);
checkloop = err<toll_d;


if  checkloop
    Psuc= Pout;
    Tsuc=Tout;
    Loopinlet=0;
    
else
    m_gas_guess= (m_gas_guess)*(1-alfa);
    counter=counter+1
end
end
end 
