function [Psuc,Tsuc] = S3PRE_inlet(p_suc,T_suc,INport_Amax,INport_Amin,V_comp,MM_g,n_van,rpm,c,toll_d,c_v,coeff_invalve,pipe,cpitch,ct,lenght,D_up,D_do,roughness,mu_g,coeff_infilter)
% This function set up the iterative approach useful for the evaluation of
%the real temperature and pressures in the inlet of the air-end section 
%or the true begining of the compression/expansion process
%
% INPUT
% P_suc [Pa]          : suction pressure at the begginning of the inlet duct
% T_suc [K]           : inlet temperature in correspondance of the inlet duct entrance
% INport_Amax [m^2]   : inlet port maximum passage area  
% INport_Amin [m^2]   : inlet port minimum passage area  
% V_comp(1) [m^3]     : first closed cell volume
% MM_g [kg/kmol]      : molecular mass of the gas 
% Rgas [J/kgK]        : ppecific gas costant 
% n_van [-]           : number of vanes
% rpm [rpm]           : shaft angular speed
% c [-]               : geometry index
%   1                 : circular stator
%   2                 : elliptical stator
% toll_d [-]          : discretization tolerance
% c_v [J/kgK]         : gas specific heat at constant volume
% aQ [kg/m^7]         : pressure loss quadratic coefficient 
% bQ [kg/m^4]         : pressure loss quadratic coefficient
% cQ [Pa]             : pressure loss quadratic coefficient
% pipe [-]            : information on type of pipe (standard,corrugated)
% cpitch [m]          : corrugated pipe pitch
% ct [m]              : corrugation height
% lenght [m]          : length of the pipe
% D_up [m]            : diameter at pipe's starting point (upstream side)
% D_do [m]            : diameter at pipe's ending point (downstream side
% roughness [micro-m] : pipe's roughness
% mu_g [Pa s]         : gas dynamic viscosity
% dQ [kg/m^7]         : pressure loss quadratic coefficient
% eQ [kg/m^4]         : pressure loss quadratic coefficient
% fQ [Pa]             : pressure loss quadratic coefficient
%
% OUTPUT
% P_suc [Pa]          : real suction pressure at the begginning of the compression process (inlet of the air-end)
% T_suc [K]           : real inlet temperature in correspondance of the air end inlet (beggining of the compression process)
% 
% HISTORY:  Gianoncelli_Genoni: creation of the file, see thesis for
% further information

 %% DEFINITIONS %%
 
 port_type="inlet";
 R_g        = SX_Constant({'UniGasConstant'})/MM_g;          % specific gas constant [J/ kg K];
 gamma= (c_v + R_g)/c_v;                                     % cp/cv
 m_gas_guess = p_suc*V_comp(1)/(R_g*T_suc)*c*n_van*rpm/60;   % guessed mass flow rate for first iteration 
 alfa=0.00001;                                               % loop variable
 Loopinlet=1;
 toll_d=1e-6;
 fb=1;
 mastseries=[0];
 while Loopinlet

[p_f,T_f,rho_f,delta_p1] = Concentrated_Losses(fb,coeff_infilter,m_gas_guess,p_suc,T_suc,MM_g,gamma);  %intake filter pressure drop

[p_f,T_f,rho_f,delta_p2] = Distributed_Losses(fb,pipe,cpitch,ct,lenght,D_up,D_do,roughness,p_f,T_f,m_gas_guess,MM_g,gamma,mu_g);

[p_f,T_f,rho_f,delta_p3] = Concentrated_Losses(fb,coeff_invalve,m_gas_guess,p_f,T_f,MM_g,gamma);  %intake valve pressure drop

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
end
 end
 
 deltap = [delta_p1 , delta_p2 , delta_p3, (p_f - Pout) , Pout]
 
end 
