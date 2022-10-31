function [p_f,T_f,rho_f,delta_p] = Distributed_Losses(fb,pipe,s,t,L,D_i,D_f,eps,p_i,T_i,m_flow,MM_g,gamma,mu_g)
% This function computes distributed pressure losses, temperature and
% density along a pipe 
% 
% INPUT
% fb      [-]           :if 1 -> forward operation, if 2 -> backward operation
% pipe    [-]           :information on type of pipe (standard,corrugated)
% s       [m]           :corrugated pipe pitch
% t       [m]           :corrugation height
% L       [m]           :length of the pipe
% D_i     [m]           :diameter at pipe's starting point
% D_f     [m]           :diameter at pipe's ending point
% eps     [micro m]     :pipe's roughness
% p_i     [Pa]          :pressure boundary condition
% T_i     [K]           :temperature boundary condition
% m_flow  [kg/s]        :mass flow rate along the pipe
% GAS                   :structure containing gas properties
% fOK                   :error flag
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
if D_i == D_f
    Re = 4*m_flow/(mu_g*pi*D_i);
    if pipe == 'standard'
       f = Darcy_TkaMil(Re,D_i,eps);
       
    elseif pipe == 'corrugated'
         if ((t/D_i) >= 0.0455) && ((t/D_i) <= 0.0635) && ((t/s) >= 0.2) &&((t/s) <= 0.6) 
           f = Darcy_Kauder(Re,D_i,t,s);
         elseif 273.15 < T_i < 313.15  
           f = Darcy_HawHelmes(D_i,s);
         else
             warning('S0_Distributed_Losses','Select a coherent model for corrugated pipe')
             SX_Logfile ('v',{lastwarn});
         end
         
    else
        warning('S0_Distributed_Losses','Select a coherent type of pipe')
        SX_Logfile ('v',{lastwarn});
       
    end
    
    delta_p = PressDrop_Incomp_Dsame(m_flow,f,D_i,L,rho_i);

else %Changing section along pipe's axial direction
    D_mean = 0.5*(D_i + D_f);                                              %Mean diameter used for Reynolds
    Re = 4*m_flow/(mu_g*pi*D_mean);
    if pipe == 'standard'
        f = Darcy_TkaMil(Re,D_mean,eps);
       
    elseif pipe == 'corrugated'
        if 0.0455 <= t/D_i <= 0.0635 && 0.2 <= t/s <= 0.6
            f = Darcy_Kauder(Re,D_mean,t,s);
        elseif 273.15 < T_i < 313.15
            f = Darcy_HawHelmes(D_mean,s);
        else
            warning('S0_Distributed_Losses','Select a coherent model for corrugated pipe')
            SX_Logfile ('v',{lastwarn});
        end
        
    else
        warning('S0_Distributed_Losses','Select a coherent type of pipe')
        SX_Logfile ('v',{lastwarn});
       
    end
    
    delta_p = PressDrop_Incomp_Ddiff(m_flow,f,D_i,D_f,L,rho_i);
    
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

function [delta_p] = PressDrop_Incomp_Ddiff(m_flow,f,D_i,D_f,L,rho)
%For circular pipes with changing diameter (linear behavior of D along
%axial direction
%Valid for incompressible fluid conditions, so Mach < 0.3
tn_Beta = tan((D_i - D_f)/(2*L));
delta_p = 1/rho * (1/(D_f^4) - 1/(D_i^4)) * (8 + f/tn_Beta)* (m_flow/pi)^2 ;
end

%% Darcy friction Coefficient Functions %%
function [f] = Darcy_TkaMil(Re,D,eps)
%For circular pipes
%Darcy friction coefficient - Tkachenko-Mileikovskyi
A0 = -0.79638*log10((eps/D)/8.208 + 7.3357/Re);
A1 = Re*eps/D + 9.3120665*A0;
f = ((8.128943 + A1)/(8.128943*A0 - 0.86859209*A1*log10(A1/3.7099535*Re)))^2;
end

function [f] = Darcy_Kauder(Re,D,t,s)
%For corrugated pipes
%Darcy friction coefficient - Kauder
%Validity conditions: ((t/D_i) >= 0.0455) && ((t/D_i) <= 0.0635) && ((t/s) >= 0.2) &&((t/s) <= 0.6)
f = 4 * exp(6.75 + 4.13 * log(t/D) + (230*((t/D)^2.1)-0.7) * log(t/s) + 0.193*exp(-3300*((t/D)^2.6)*(t/s)) * log(Re)); 
end

function [f] = Darcy_HawHelmes(D,s)
%For corrugated pipes
%Darcy friction coefficient - Hawthorne and Von Helmes 
%Valid for air and water at ambient temperature 
f = (D/s)*(1-(D/(D + 0.438*s))^2)^2;       
end



        
    
    