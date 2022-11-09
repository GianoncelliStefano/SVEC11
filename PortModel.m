function [ Pout , Tout] = PortModel(Pin , Tin , Amax , Amin, gamma, Rgas,  m, port_type)
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
