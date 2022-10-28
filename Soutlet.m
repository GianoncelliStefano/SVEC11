function [ Pout , Tout] = Soutlet( Pin , Tin , OUTport_Amax, OUTport_Amin,mgas, MM_g,c_v, Pgeom )
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
 port_type="outlet";
 R_g        = SX_Constant({'UniGasConstant'})/MM_g;  % specific gas constant [J/ kg K];
 gamma= (c_v + R_g)/c_v;
 alfa=0.001;
 Pd=Pin; %Sara poi da aggiornare aggiungendo la possibilita di mettere la pressione di genny
 LoopOutlet=1;
 toll_d=1e-6;
 while LoopOutlet
[m_port_inf , Uinf , Pinf ,Tinf, Pout , Tout] = PortModel(Pd , Tin , OUTport_Amax , OUTport_Amin, gamma, R_g,  mgas, port_type);
err= abs(Pout-Pgeom);
checkloop = err<toll_d;
if  checkloop
    Loopinlet=0;
else
    Pd = Pd*(1+alfa);
end
end
end 
