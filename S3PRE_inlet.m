function [p_in,T_in,deltap_inlet] = S3PRE_inlet(p_suc,T_suc,INport_Amax,INport_Amin,V_comp,MM_g,n_van,rpm,c,c_v,coeff_invalve,pipe,cpitch,ct,lenght,D_up,D_do,roughness,mu_g,coeff_infilter,fSDP)
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
% fSDP                : suction and discharge model activation flag
%
% OUTPUT
% P_suc [Pa]          : real suction pressure at the begginning of the compression process (inlet of the air-end)
% T_suc [K]           : real inlet temperature in correspondance of the air end inlet (beggining of the compression process)
% 
% HISTORY:  Gianoncelli_Genoni: creation of the file, see thesis for
% further information

 %% fSDP flag check %%
 switch fSDP
     case 0 
         p_in         = p_suc;                                % when the intake model is not active, the air-end inlet pressure is equal to the suction one
         T_in         = T_suc;                                % when the intake model is not active, the air-end inlet temperature is equal to the suction one
         deltap_inlet = 0;                                    % when the intake model is not active, the pressure drop during the suction process is null 
     case 1
 port_type   = "inlet";
 R_g         = SX_Constant({'UniGasConstant'})/MM_g;          % specific gas constant [J/ kg K];
 gamma       = (c_v + R_g)/c_v;                               % heat capacity ratio
 m_gas_guess = p_suc*V_comp(1)/(R_g*T_suc)*c*n_van*rpm/60;    % guessed mass flow rate for first iteration 
 alfa        = 0.00001;                                       % under-relaxation factor for the inlet loop
 Loopinlet   = 1;                                             % loop control variable 
 toll_d      = 1e-6;                                          % mass flow rate tolerance for the loop
 fb          = 1;                                             % flag for forward computation of pressure losses
 
 while Loopinlet

[p_df,T_df,delta_pfilter] = Concentrated_Losses(fb,coeff_infilter,m_gas_guess,p_suc,T_suc,MM_g,gamma);                               % intake filter concentrated pressure drop

[p_dp,T_dp,delta_pduct]   = Distributed_Losses(fb,pipe,cpitch,ct,lenght,D_up,D_do,roughness,p_df,T_df,m_gas_guess,MM_g,gamma,mu_g);  % intake duct distributed pressure drop

[p_dv,T_dv,delta_pvalve]  = Concentrated_Losses(fb,coeff_invalve,m_gas_guess,p_dp,T_dp,MM_g,gamma);                                  % intake valve concentrated pressure drop

[Pout , Tout]             = PortModel(p_dv , T_dv , INport_Amax , INport_Amin, gamma, R_g,  m_gas_guess, port_type);                 % inlet port pressure drop

mast      = Pout*V_comp(1)/(R_g*Tout)*c*n_van*rpm/60;     % updated mass flow rate                                                                           
err       = abs(mast-m_gas_guess);                        % residual error on the control variable
checkloop = err < toll_d;                                 % loop exit condition 


if  checkloop
    p_in      = Pout;                                     % when converged, the air-end inlet pressure is equal to the loop output one
    T_in      = Tout;                                     % when converged, the air-end inlet temperature is equal to the loop output one
    Loopinlet = 0;                                          
    
else
    m_gas_guess = (m_gas_guess)*(1-alfa);                 % loop decisional variable update through under-relaxation factor
end
 end
 
 deltap_inlet = delta_pfilter + delta_pduct + delta_pvalve +  (p_dv - Pout);     % intake process overall pressure loss 
 deltap       = [delta_pfilter , delta_pduct , delta_pvalve, p_dv - Pout, Pout]  % intake process partial pressure loss on each component

 end

 %% INTERNAL FUNCTIONS  %%

function [p_f,T_f,delta_p] = Concentrated_Losses(fb,coeff,m_flow,p_i,T_i,MM_g,gamma)
% This function computes concentrated pressure losses, temperature and
% density across a flow disturbing element
%
% INPUT
% fb      [-]          :if 1 -> forward operation, if 2 -> backward operation
% coeff   [-]          :structure containing pressure loss quadratic coefficients
% m_flow  [kg/s]       :mass flow rate at the element
% p_i     [Pa]         :pressure boundary condition
% T_i     [K]          :temperature boundary condition
% MM_g    [kg/mol]     :molar mass of gas 
% gamma   [-]          :cp/cv
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
% else
%    warning('S0_Concentrated_Losses','Wrong number inserted for fb, please insert 1 for forward computing or 2 for backward computing')
%    SX_Logfile ('v',{lastwarn});
end

%% Concentrated Losses Computing %%
delta_p = coeff(1)*(Q^2) + coeff(2)*Q + coeff(3);
p_f = p_i - delta_p*coeffDP;
T_f = T_i * (p_f/p_i)^((gamma-1)/gamma);            %HP: adiabatic element
rho_f = p_f*MM_g/(R_u*T_f);

end
function [p_f,T_f,delta_p] = Distributed_Losses(fb,pipe,s,t,L,D_up,D_do,eps,p_i,T_i,m_flow,MM_g,gamma,mu_g)

% This function computes distributed pressure losses, temperature and
% density along a pipe 
% 
% INPUT
% fb      [-]           :if 1 -> forward operation, if 2 -> backward operation
% pipe    [-]           :information on type of pipe (standard,corrugated)
% s       [m]           :corrugated pipe pitch
% t       [m]           :corrugation height
% L       [m]           :length of the pipe
% D_up    [m]           :diameter at pipe's starting point (upstream side)
% D_do    [m]           :diameter at pipe's ending point (downstream side)
% eps     [micro m]     :pipe's roughness
% p_i     [Pa]          :pressure boundary condition
% T_i     [K]           :temperature boundary condition
% m_flow  [kg/s]        :mass flow rate along the pipe
% MM_g    [kg/mol]      :molar mass of gas 
% gamma   [-]           :cp/cv
% mu_g    [Pa s]        :gas dynamic viscosity
%
% OUTPUT
% p_f     [Pa]          :resulting pressure
% T_f     [K]           :resulting temperature
% rho_f   [kg/m^3]      :resulting density
% delta_p [Pa]          :resulting pressure drop
%
% Developers: Genoni-Gianoncelli

%% Starting Conditions %%
R_u    = 8.314472;           
rho_i  = p_i*MM_g/(R_u*T_i);     

%% Forward or Backward Process %%
if fb == 1                         %forward -> pressure loss | backward -> pressure gain 
   coeffDP = 1;
elseif fb == 2
   coeffDP = -1;
% else
%    warning('S0_Distributed_Losses','Wrong number inserted for fb, please insert 1 for forward computing or 2 for backward computing')
%    SX_Logfile ('v',{lastwarn});
end

%% Distributed losses and friction coefficient computing %%
if L == 0
    p_f     = p_i;
    T_f     = T_i;
    rho_f   = rho_i;
    delta_p = 0;
elseif D_up == D_do
    Re = 4*m_flow/(mu_g*pi*D_up);
    if pipe == 'standard'
       f = Darcy_SwameeJain(Re,D_up,eps);
       
    elseif pipe == 'corrugated'
         if ((t/D_up) >= 0.0455) && ((t/D_up) <= 0.0635) && ((t/s) >= 0.2) &&((t/s) <= 0.6) 
           f = Darcy_Kauder(Re,D_up,t,s);
         elseif 273.15 < T_i < 313.15  
           f = Darcy_HawHelmes(D_up,s,t);
         else
             warning('S0_Distributed_Losses','Select a coherent model for corrugated pipe')
             SX_Logfile ('v',{lastwarn});
         end
         
    else
        warning('S0_Distributed_Losses','Select a coherent type of pipe')
        SX_Logfile ('v',{lastwarn});
       
    end
    
    delta_p = PressDrop_Incomp_Dsame(m_flow,f,D_up,L,rho_i);

else %Changing section along pipe's axial direction
    D_mean = 0.5*(D_up + D_do);                                              %Mean diameter used for Reynolds
    Re = 4*m_flow/(mu_g*pi*D_mean);
    if pipe == 'standard'
        f = Darcy_TkaMil(Re,D_mean,eps);
       
    elseif pipe == 'corrugated'
        if 0.0455 <= t/D_up <= 0.0635 && 0.2 <= t/s <= 0.6
            f = Darcy_Kauder(Re,D_mean,t,s);
        elseif 273.15 < T_i < 313.15
            f = Darcy_HawHelmes(D_mean,s,t);
        else
            warning('S0_Distributed_Losses','Select a coherent model for corrugated pipe')
            SX_Logfile ('v',{lastwarn});
        end
        
    else
        warning('S0_Distributed_Losses','Select a coherent type of pipe')
        SX_Logfile ('v',{lastwarn});
       
    end
    
    delta_p = PressDrop_Incomp_Ddiff(m_flow,f,D_up,D_do,L,rho_i);
    
end

p_f = p_i - delta_p*coeffDP;
T_f = T_i * (p_f/p_i)^((gamma-1)/gamma);                                   %HP: adiabatic pipe
rho_f = p_f*MM_g/(R_u*T_f);

%% Pressure Drop Functions %%
function [delta_p] = PressDrop_Incomp_Dsame(m_flow,f,D,L,rho)
%For circular pipes with constant diameter
%Valid for incompressible fluid conditions, so Mach < 0.3
delta_p = 8*f*(L/rho) * (1/(D^5)) * (m_flow/pi)^2;
end

function [delta_p] = PressDrop_Incomp_Ddiff(m_flow,f,D_up,D_do,L,rho)
%For circular pipes with changing diameter (linear behavior of D along
%axial direction
%Valid for incompressible fluid conditions, so Mach < 0.3
tn_Beta = tan((D_up - D_do)/(2*L));
delta_p = 1/rho * (1/(D_do^4) - 1/(D_up^4)) * (8 + f/tn_Beta)* (m_flow/pi)^2 ;
end

%% Darcy friction Coefficient Functions %%
function [f] = Darcy_TkaMil(Re,D,eps)
%For circular pipes
%Darcy friction coefficient - Tkachenko-Mileikovskyi
A0 = -0.79638*log10((eps/D)/8.208 + 7.3357/Re);
A1 = Re*eps/D + 9.3120665*A0;
f = ((8.128943 + A1)/(8.128943*A0 - 0.86859209*A1*log10(A1/3.7099535*Re)))^2;
end

function [f] = Darcy_SwameeJain(Re,D,eps)
%For circular pipes 
%Darcy friction coefficient - Swamee-Jain 
f = 0.25/((log10((eps/(3.7*D))+5.74/(Re^0.9)))^2);
end

function [f] = Darcy_Kauder(Re,D,t,s)
%For corrugated pipes
%Darcy friction coefficient - Kauder
%Validity conditions: ((t/D_up) >= 0.0455) && ((t/D_up) <= 0.0635) && ((t/s) >= 0.2) &&((t/s) <= 0.6)
f = 4 * exp(6.75 + 4.13 * log(t/D) + (230*((t/D)^2.1)-0.7) * log(t/s) + 0.193*exp(-3300*((t/D)^2.6)*(t/s)) * log(Re)); 
end

function [f] = Darcy_HawHelmes(D,s,t)
%For corrugated pipes
%Darcy friction coefficient - Hawthorne and Von Helmes 
%Valid for air and water at ambient temperature
d = D - 2*t; %internal diameter
f = (d/s)*(1-(d/(d + 0.438*s))^2)^2;       
end
end
function [ Pout , Tout]    = PortModel(Pin,Tin,Amax,Amin,gamma,Rgas,m,port_type)
% This function evaluate with a simplified approach the gas-dynamic 
% effects happening in correspondance of the inlet and outlet ports.
% The mean fluid velocity at steady state will be computed along with 
% the steady state pressure in the middle of the port itself
%
% INPUT
% Pin: Input pressure[Pa] to the port model i.e. the known one which is:
%       - The UPSTREAM one in the inlet port case
%       - The DOWNSTREAM one in the outlet port case
% Tin: Input temperature [K] to the port model  i.e. the known one which is: 
%      - The UPSTREAM one in the inlet port case
%      - The DOWNSTREAM one in the outlet port case
% Amax: Maximum port area [m^2]
% Amin: Minimum port area [m^2]
% gamma: Specific heat ratio of the working fluid [-]
% Rgas: Specific gas costant of the working fluid [J/kgK]
% m: mass flow rate first guess (case inlet port) or mass flow rate value(in the case of outlet port) 
% port_type: Flag useful to access the correct model depending on the port
% type (inlet or outlet)
% 
% 
% 
% 
% OUTPUT
% m_port_inf: fluid mass flow rate passing through the port once steady state is reached [kg/s]
% Uinf: Velocity of the fluid in the middle point of the port once steady state is reached [m/s]
% Pinf: Pressure of the fluid in the port middle point once steady state is reached [Pa]
% Tinf: Temperature of the fluid in the port middle point once steady state is reached [Pa]
% Pout: Output pressure [Pa] to the port model i.e the unknown one which is: 
%       - The UPSTREAM one in the inlet port case
%       - The DOWNSTREAM one in the outlet port case
% Tout: Output temperature [Pa] from the port model i.e the unknown one which is: 
%       - The UPSTREAM one in the inlet port case
%       - The DOWNSTREAM one in the outlet port case
% 
% 
% NOTE : Model based on the paper "Theoretical modeling and experimental investigations for the improvement
% of the mechanical efficiency in sliding vane rotary compressors" by G. Bianchi and R. Cipollone
% 
% HISTORY:  Gianoncelli_Genoni: creation of the file

%% DEFINITIONS %%

Xi = Amin/Amax;                                 %Contraction ratio [-]
Amean=mean([Amax,Amin]);                        %average port area [m^2]

switch port_type
    case "inlet"
        [Pinf,Tinf, Uinf,Pout, Tout]=Inlet_simplemodel(Pin, Tin, Xi,Amean,gamma,Rgas,m);
    case "outlet"
        [Pinf, Tinf,Uinf,Pout,Tout]=Outlet_simplemodel(Pin,Tin,Xi,Amean,gamma,Rgas,m);
end

m_port_inf = Pinf*Uinf*Amean/(Rgas*Tinf);       %steady state mass flow rate [kg/s]  



%% MODEL SOLUTION AT STEADY STATE %%
function [Pss, Tss, Uss,Pd, Td]=Inlet_simplemodel(P,T,Xi,Amean,gamma,R,m)
fun1= @(Pd) P/(R*T)*(1-Xi^2*(1-(Pd/P)^((gamma-1)/gamma)))^(1/(gamma-1))*Xi*Amean*(((2*gamma*R*T)/(gamma-1))*(1-(Pd/P)^((gamma-1)/gamma)))^(1/2)-m;
Pdrange=[10000 P];
Pd=fzero(fun1,Pdrange);
Uss=((2*gamma*R*T)/(gamma-1)*Xi^2*(1-(Pd/P)^((gamma-1)/gamma)))^(1/2);
Pss=P*(1-Xi^2*(1-(Pd/P)^((gamma-1)/gamma)))^(gamma/(gamma-1));
Tss=T*(Pss/P)^((gamma-1)/gamma);
Td=T*(Pd/P)^((gamma-1)/gamma);
end

function [Pss , Tss, Uss, Pu , Tu]=Outlet_simplemodel(P,T,Xi, Amean,gamma,R,m)
Tu=T;
fun2= @(Pu) Pu/(R*Tu)*(1-Xi^2*(1-(P/Pu)^((gamma-1)/gamma)))^(1/(gamma-1))*Xi*Amean*(((2*gamma*R*Tu)/(gamma-1))*(1-(P/Pu)^((gamma-1)/gamma)))^(1/2)-m;
Purange=[P 50e5];
Pu=fzero(fun2,Purange);
Uss=((2*gamma*R*Tu)/(gamma-1)*Xi^2*(1-(P/Pu)^((gamma-1)/gamma)))^(1/2);
Pss=Pu*(1-Xi^2*(1-(P/Pu)^((gamma-1)/gamma)))^(gamma/(gamma-1));
Tss=Tu*(Pss/Pu)^((gamma-1)/gamma);
end
end 

end 
