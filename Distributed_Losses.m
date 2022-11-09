function [p_f,T_f,rho_f,delta_p] = Distributed_Losses(fb,pipe,s,t,L,D_up,D_do,eps,p_i,T_i,m_flow,MM_g,gamma,mu_g)
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
eps    = eps * 10^-6;
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

end

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



        
    
    