function [p_in,T_in,deltap_inlet] = S3_PRE_inlet(p_suc,T_suc,INport_Amax,INport_Amin,V_comp1,MM_g,n_van,rpm,c,c_v,coeff_invalve,pipe,cpitch,ct,lenght,D_up,D_do,roughness,mu_g,coeff_infilter,fSDP)
% This function sets up the iterative approach useful for the evaluation of
% the actual temperature and pressure at the inlet of the air-end section 
%
% INPUT
% p_suc          [Pa]       : ambient pressure
% T_suc          [K]        : ambient temperature 
% INport_Amax    [m^2]      : inlet port maximum passage area  
% INport_Amin    [m^2]      : inlet port minimum passage area  
% V_comp1        [m^3]      : first closed cell volume
% MM_g           [kg/kmol]  : molecular mass of the gas  
% n_van          [-]        : number of vanes
% rpm            [rpm]      : shaft angular speed
% c              [-]        : geometry index
%   1                       : circular stator
%   2                       : elliptical stator
% c_v            [J/kgK]    : gas specific heat at constant volume
% coeff_invalve  [-]        : intake valve pressure losses coefficients
% pipe           [-]        : information on type of pipe (standard,corrugated)
% cpitch         [m]        : corrugated pipe pitch
% ct             [m]        : corrugation height
% lenght         [m]        : length of the pipe
% D_up           [m]        : diameter at pipe's starting point (upstream side)
% D_do           [m]        : diameter at pipe's ending point (downstream side)
% roughness      [micro-m]  : pipe's roughness
% mu_g           [Pa s]     : gas dynamic viscosity
% coeff_infilter [-]        : filter pressure losses coefficients
% fSDP           [-]        : suction and discharge model activation flag
%
% OUTPUT
% p_in           [Pa]       : inlet pressure at the begginning of the closed cell process (inlet of the air-end)
% T_in           [K]        : inlet temperature at the begginning of the closed cell process (inlet of the air-end)
% deltap_inlet   [Pa]       : inlet process pressure loss (p_in - p_suc)
% 
% HISTORY:  Gianoncelli_Genoni: creation of the file, see thesis for further information

 %% INTAKE PROCESS LOOP %%
 switch fSDP
    case 0                                                   % intake process model not activated
         p_in         = p_suc;                                % when the intake model is not active, the air-end inlet pressure is equal to the suction one
         T_in         = T_suc;                                % when the intake model is not active, the air-end inlet temperature is equal to the suction one
         deltap_inlet = 0;                                    % when the intake model is not active, the pressure drop during the suction process is null 
    case 1                                                   % intake process model activated
 port_type   = "inlet";                                       
 R_g         = SX_Constant({'UniGasConstant'})/MM_g;          % specific gas constant [J/kg K];
 gamma       = (c_v + R_g)/c_v;                               % heat capacity ratio
 m_gas_guess = p_suc*V_comp1/(R_g*T_suc)*c*n_van*rpm/60;      % guessed mass flow rate for first iteration 
 alfa        = 1e-5;                                          % under-relaxation factor for the inlet loop
 Loopinlet   = 1;                                             % loop control variable 
 toll_d      = 1e-6;                                          % mass flow rate tolerance for the loop
 fb          = 1;                                             % flag for forward computation of pressure losses
 
  while Loopinlet

[p_df,T_df,delta_pfilter]  = Concentrated_Losses(fb,coeff_infilter,m_gas_guess,p_suc,T_suc,MM_g,gamma);                               % intake filter concentrated pressure drop

[p_dp,T_dp,delta_pduct]    = Distributed_Losses(fb,pipe,cpitch,ct,lenght,D_up,D_do,roughness,p_df,T_df,m_gas_guess,MM_g,gamma,mu_g);  % intake duct distributed pressure drop

[p_dv,T_dv,delta_pvalve]   = Concentrated_Losses(fb,coeff_invalve,m_gas_guess,p_dp,T_dp,MM_g,gamma);                                  % intake valve concentrated pressure drop

[Pout , Tout, delta_pport] = PortModel(p_dv , T_dv , INport_Amax , INport_Amin, gamma, R_g,  m_gas_guess, port_type);                 % inlet port pressure drop

mast      = Pout*V_comp1/(R_g*Tout)*c*n_van*rpm/60;     % updated mass flow rate                                                                           
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
 
 deltap_inlet = delta_pfilter + delta_pduct + delta_pvalve +  delta_pport;     % intake process overall pressure loss 
 deltap       = [delta_pfilter , delta_pduct , delta_pvalve, delta_pport , Pout]  % intake process partial pressure loss on each component

 end
end

%% INTERNAL FUNCTIONS  %%

function [p_f,T_f,delta_p]   = Concentrated_Losses(fb,coeff,m_flow,p_i,T_i,MM_g,gamma)
% This function computes concentrated pressure losses and temperature variation across a flow disturbing element
%
% INPUT
% fb      [-]          :if 1 -> forward operation, if 2 -> backward operation
% coeff   [-]          :structure containing pressure loss quadratic coefficients
% m_flow  [kg/s]       :mass flow rate at the element
% p_i     [Pa]         :pressure boundary condition
% T_i     [K]          :temperature boundary condition
% MM_g    [kg/mol]     :molar mass of gas 
% gamma   [-]          :heat capacity ratio
%
% OUTPUT
% 
% p_f     [Pa]         :resulting pressure
% T_f     [K]          :resulting temperature
% delta_p [Pa]         :resulting pressure drop
%
% Developers: Genoni-Gianoncelli

% Starting Conditions %
R_u    = SX_Constant({'UniGasConstant'});           % universal gas constant [J/mol K]
rho_i  = p_i*MM_g/(R_u*T_i);                        % density at element beginning [kg/m3]
Q      = m_flow/rho_i;                              % volumetric flow rate [m3/s]

% Forward or Backward Process %
switch fb
    case 1                                          %forward --> pressure loss | backward --> pressure gain 
   coeffDP = 1;
    case 2
   coeffDP = -1;
end

% Concentrated Losses Computation %
delta_p = coeff(1)*(Q^2) + coeff(2)*Q + coeff(3);
p_f = p_i - delta_p*coeffDP;
T_f = T_i * (p_f/p_i)^((gamma-1)/gamma);            %HP: adiabatic element
%rho_f = p_f*MM_g/(R_u*T_f);

 end
function [p_f,T_f,delta_p]   = Distributed_Losses(fb,pipe,s,t,L,D_up,D_do,eps,p_i,T_i,m_flow,MM_g,gamma,mu_g)
% This function computes returns pressure, temperature and distributed
% pressure losses along a pipe
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
% delta_p [Pa]          :resulting pressure drop
%
% Developers: Genoni-Gianoncelli

% checks %
check1 = L ~= 0;                                         % check if an intake pipe is present: 1 --> yes 0 --> no
check2 = L>0 && D_up==D_do;                              % check if the section is costant or not 1 --> costant section 0 --> variable section

switch check1

case 0                                                   % no intake pipe present
    p_f     = p_i;
    T_f     = T_i;
    delta_p = 0;
case 1                                                   % intake pipe is present
% Starting Conditions %
R_u    = SX_Constant({'UniGasConstant'});                % universal gas constant [J/mol K] 
rho_i  = p_i*MM_g/(R_u*T_i);                             % density at pipe beginning [kg/m3]

% Forward or Backward Process %
  switch fb
        case 1                                           %forward -> pressure loss | backward -> pressure gain 
             coeffDP = 1;
        case 2
             coeffDP = -1;
  end

% Distributed losses and friction coefficient computing %
switch check2
case 1                                                   % costant section                           
            Re = 4*m_flow/(mu_g*pi*D_up);                
         if pipe == "standard"
            f = Darcy_SwameeJain(Re,D_up,eps);
         elseif pipe == "corrugated"
                if ((t/D_up) >= 0.0455) && ((t/D_up) <= 0.0635) && ((t/s) >= 0.2) &&((t/s) <= 0.6) 
                    f = Darcy_Kauder(Re,D_up,t,s);
                elseif 273.15 < T_i < 313.15  
                    f = Darcy_HawHelmes(D_up,s);
                else
             warning('S3PRE_inlet: Select a coherent model for corrugated pipe')
             SX_Logfile ('v',{lastwarn});
                end
         
          else
        warning('S3PRE_inlet: Select a coherent type of pipe')
        SX_Logfile ('v',{lastwarn});
       
          end
    
delta_p = PressDrop_Incomp_Dsame(m_flow,f,D_up,L,rho_i);   % pressure drop computation for ducts with costant section

case 0                                                     % variable section 
    D_mean = 0.5*(D_up + D_do);                          
    Re = 4*m_flow/(mu_g*pi*D_mean);
    if pipe == "standard"
        f = Darcy_SwameeJain(Re,D_mean,eps);
       
    elseif pipe == "corrugated"
        if ((t/D_up) >= 0.0455) && ((t/D_up) <= 0.0635) && ((t/s) >= 0.2) &&((t/s) <= 0.6) 
            f = Darcy_Kauder(Re,D_mean,t,s);
        elseif (273.15 < T_i) && (T_i < 313.15)
            f = Darcy_HawHelmes(D_mean,s);
        else
            warning('S3PRE_inlet: Select a coherent model for corrugated pipe')
            SX_Logfile ('v',{lastwarn});
        end
        
    else
        warning('S3PRE_inlet: Select a coherent type of pipe')
        SX_Logfile ('v',{lastwarn});
       
    end
    
delta_p = PressDrop_Incomp_Ddiff(m_flow,f,D_up,D_do,L,rho_i);% pressure drop computation for ducts with variable section
    
end
end

p_f = p_i - delta_p*coeffDP;
T_f = T_i * (p_f/p_i)^((gamma-1)/gamma);                     %HP: adiabatic pipe
% rho_f = p_f*MM_g/(R_u*T_f);

%% Pressure Drop Functions %%
function [delta_p] = PressDrop_Incomp_Dsame(m_flow,f,D,L,rho)
%For circular pipes with constant diameter
%Valid for incompressible fluid conditions (Mach < 0.3)
delta_p = 8*f*(L/rho) * (1/(D^5)) * (m_flow/pi)^2;
end

function [delta_p] = PressDrop_Incomp_Ddiff(m_flow,f,D_up,D_do,L,rho)
%For circular pipes with changing diameter (linear behavior of D along
%axial direction
%Valid for incompressible fluid conditions (Mach < 0.3)
tn_Beta = ((D_up - D_do)/(2*L));
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

function [f] = Darcy_HawHelmes(D,s)
%For corrugated pipes
%Darcy friction coefficient - Hawthorne and Von Helmes 
%Valid for air and water at ambient temperature
f = (D/s)*(1-(D/(D + 0.438*s))^2)^2;       
end

end
function [Pout,Tout,delta_p] = PortModel(Pin,Tin,Amax,Amin,gamma,Rgas,m,port_type)
% This function evaluates, with a simplified approach, the gas-dynamic 
% effects happening in correspondance of the inlet and outlet ports.
% The mean fluid velocity at steady state will be computed along with 
% the steady state pressure in the middle of the port itself
%
% INPUT
% Pin        [Pa]     : Input pressure to the port model i.e. the known one which is:
%                       - The UPSTREAM one in the inlet port case
%                       - The DOWNSTREAM one in the outlet port case
% Tin        [K]      : Input temperature  to the port model  i.e. the known one which is: 
%                       - The UPSTREAM one in the inlet port case
%                       - The DOWNSTREAM one in the outlet port case
% Amax       [m^2]    : Maximum port area 
% Amin       [m^2]    : Minimum port area 
% gamma      [-]      : Specific heat ratio 
% Rgas       [J/kgK]  : Specific gas costant 
% m          [kg/s]   : mass flow rate first guess (case inlet port) or mass flow rate value (in the case of outlet port) 
% port_type  [-]      : Flag useful to access the correct model depending on the port type (inlet or outlet) 
% 
% OUTPUT
% Pout       [Pa]     : Output pressure to the port model i.e the unknown one which is: 
%                       - The UPSTREAM one in the inlet port case
%                       - The DOWNSTREAM one in the outlet port case
% Tout       [Pa]     : Output temperature from the port model i.e the unknown one which is: 
%                       - The UPSTREAM one in the inlet port case
%                       - The DOWNSTREAM one in the outlet port case
% 
% ADDITIONAL
% m_port_inf [kg/s]   : fluid mass flow rate passing through the port once steady state is reached 
% Uinf       [m/s]    : Velocity of the fluid in the middle point of the port once steady state is reached 
% Pinf       [Pa]     : Pressure of the fluid in the port middle point once steady state is reached 
% Tinf       [Pa]     : Temperature of the fluid in the port middle point once steady state is reached 
% 
% NOTE : Model based on the paper "Theoretical modeling and experimental investigations for the improvement
% of the mechanical efficiency in sliding vane rotary compressors" by G. Bianchi and R. Cipollone
% 
% HISTORY:  Gianoncelli_Genoni: creation of the file

% DEFINITIONS %
Xi = Amin/Amax;                                 % contraction ratio [-]
Amean=mean([Amax,Amin]);                        % average port area [m^2]

switch port_type
    case "inlet"
        [Pinf,Tinf, Uinf,Pout, Tout] = Inlet_PortModel(Pin, Tin, Xi,Amean,gamma,Rgas,m);
    case "outlet"
        [Pinf, Tinf,Uinf,Pout,Tout]  = Outlet_PortModel(Pin,Tin,Xi,Amean,gamma,Rgas,m);
end

m_port_inf = Pinf*Uinf*Amean/(Rgas*Tinf);       % steady state mass flow rate [kg/s]  
delta_p = Pin - Pout;                             % pressure drop over the port

% MODEL SOLUTION AT STEADY STATE %
function [Pss,Tss,Uss,Pd,Td] = Inlet_PortModel(P,T,Xi,Amean,gamma,R,m)
fun1= @(Pd) P/(R*T)*(1-Xi^2*(1-(Pd/P)^((gamma-1)/gamma)))^(1/(gamma-1))*Xi*Amean*(((2*gamma*R*T)/(gamma-1))*(1-(Pd/P)^((gamma-1)/gamma)))^(1/2)-m;
Pdrange=[10000 P];
Pd=fzero(fun1,Pdrange);
Uss=((2*gamma*R*T)/(gamma-1)*Xi^2*(1-(Pd/P)^((gamma-1)/gamma)))^(1/2);
Pss=P*(1-Xi^2*(1-(Pd/P)^((gamma-1)/gamma)))^(gamma/(gamma-1));
Tss=T*(Pss/P)^((gamma-1)/gamma);
Td=T*(Pd/P)^((gamma-1)/gamma);
end

function [Pss,Tss,Uss,Pu,Tu] = Outlet_PortModel(P,T,Xi, Amean,gamma,R,m)
Tu=T;
fun2= @(Pu) Pu/(R*Tu)*(1-Xi^2*(1-(P/Pu)^((gamma-1)/gamma)))^(1/(gamma-1))*Xi*Amean*(((2*gamma*R*Tu)/(gamma-1))*(1-(P/Pu)^((gamma-1)/gamma)))^(1/2)-m;
Purange=[P 50e5];
Pu=fzero(fun2,Purange);
Uss=((2*gamma*R*Tu)/(gamma-1)*Xi^2*(1-(P/Pu)^((gamma-1)/gamma)))^(1/2);
Pss=Pu*(1-Xi^2*(1-(P/Pu)^((gamma-1)/gamma)))^(gamma/(gamma-1));
Tss=Tu*(Pss/Pu)^((gamma-1)/gamma);
end

end 
 
