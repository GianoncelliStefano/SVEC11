function [p_f,T_f,rho_f,delta_p] = Concentrated_Losses(fb,coeff,m_flow,p_i,T_i,MM_g,gamma)
% This function computes concentrated pressure losses, temperature and
% density across a flow disturbing element
%
% INPUT
% fb      [-]          :if 1 -> forward operation, if 2 -> backward operation
% ELEM    [-]          :structure containing pressure loss quadratic coefficients
% m_flow  [kg/s]       :mass flow rate at the element
% p_i     [Pa]         :pressure boundary condition
% T_i     [K]          :temperature boundary condition
% D       [m]          :diameter of the element
% GAS                  :structure containing gas properties
% fOK                  :error flag          
%
% OUTPUT
% 
% p_f     [Pa]         :resulting pressure
% T_f     [K]          :resulting temperature
% rho_f   [kg/m^3]     :resulting density
% delta_p [Pa]         :resulting pressure drop
%
% Developers: Genoni-Gianoncelli

%% Starting Conditions %%
R_u    = 8.314472;                 
rho_i  = p_i*MM_g/(R_u*T_i);     
Q      = m_flow/rho_i;                              %volumetric flow rate [m3/s]

%% Forward or Backward Process %%
if fb == 1                                          %forward --> pressure loss | backward --> pressure gain 
   coeffDP = 1;
elseif fb == 2
   coeffDP = -1;
else
   warning('S0_Concentrated_Losses','Wrong number inserted for fb, please insert 1 for forward computing or 2 for backward computing')
   SX_Logfile ('v',{lastwarn});
end

%% Concentrated Losses Computing %%
delta_p = coeff(1)*(Q^2) + coeff(2)*Q + coeff(3);
p_f = p_i - delta_p*coeffDP;
T_f = T_i * (p_f/p_i)^((gamma-1)/gamma);            %HP: adiabatic element
rho_f = p_f*MM_g/(R_u*T_f);

end

