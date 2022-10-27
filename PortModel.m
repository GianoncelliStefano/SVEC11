function [m_port_inf , Uinf , Pinf ,Tinf, Pout , Tout] = PortModel(Pin , Tin , Amax , Amin, gamma, Rgas,  m, port_type)
% This function evaluate with a simplified approach the gas-dynamic 
% effects happening in correspondance of the inlet and outlet ports
% The mean fluid velocity at steady state will be computed along with 
% the average steady state pressure in the port itself
%
% INPUT
% Pin: Input pressure[Pa] to the port model i.e. the known one which is:
%       - The UPSTREAM one in the inlet port case
%       - The DOWNSTREAM one in the outlet port case
% Tin: Input temperature [K] to the port model  i.e. the known one which is: 
%      - The UPSTREAM one in the inlet port case
%      - The DOWNSTREAM one in the outlet port case
% Au: Upstream port area [m^2]
% Ad: DOwnstream port area [m^2]
% gamma: Specific heat ratio of the fluid passing through the port [-]
% Rgas: Specific gas costant [J/kgK]
% m: mass flow rate first guess (case inlet port) or mass flow rate value(in the case of outlet port) 
% port_type: Flag useful to access the correct model depending on the port
% under investigation: the inlet or the outlet one
% Tsvec: Temperature at the end of the compression (closed cell), assumed
% equal to the one in the first part of the outlet duct
% 
% 
% 
% OUTPUT
% m_port_inf: fluid mass flow rate passing through the port once steady state is reached [kg/s]
% Uinf: Velocity of the fluid in the middle point of the port once steady state is reached [m/s]
% Pinf: Pressure of the fluid in the port middle point once steady state is reached [Pa]
% Tinf: Temperature of the fluid in the port middle point once steady state is reached [Pa]
% Pout: Output pressure [Pa] to the port model i.e the unknown one which is: 
%       -The UPSTREAM one in the inlet port case
%       - The DOWNSTREAM one in the outlet port case
% Tout: Output temperature [Pa] from the port model i.e the unknown one which is: 
%       - The UPSTREAM one in the inlet port case
%       - The DOWNSTREAM one in the outlet port case
% 
% 
% NOTE : Model based on the paper "Theoretical modeling and experimental investigations for the improvement
% of the mechanical efficiency in sliding vane rotary compressors" by G. Bianchi and R. Cipollone
% 
% HISTORY:  Gianoncelli_Genoni: creation of the file and addition of the port model in SVEC

    %% DEFINITIONS %%

Xi = Amin/Amax; %Contraction ratio [-]
Amean=mean([Amax,Amin]); %average port area [m^2]

switch port_type
    case "inlet"
        [Pinf,Tinf, Uinf,Pout, Tout]=Inlet_simplemodel(Pin, Tin, Xi,Amean,gamma,Rgas,m);
    case "outlet"
        [Pinf, Tinf,Uinf,Pout,Tout]=Outlet_simplemodel(Pin,Tin,Xi,Amean,gamma,Rgas,m);UGGYTYT
end
m_port_inf = Pinf*Uinf*Amean/(Rgas*Tinf);



%% MODEL SOLUTION AT STEADY STATE - simple model %%
function [Pss, Tss, Uss,Pd, Td]=Inlet_simplemodel(P,T,Xi,Amean,gamma,R,m)
fun1= @(Pd) P/(R*T)*(1-Xi^2*(1-(Pd/P)^((gamma-1)/gamma)))^(1/(gamma-1))*Xi*Amean*(((2*gamma*R*T)/(gamma-1))*(1-(Pd/P)^((gamma-1)/gamma)))^(1/2)-m;
Pdrange=[0.1 P];
Pd=fzero(fun1,Pdrange);
Uss=((2*gamma*R*T)/(gamma-1)*Xi^2*(1-(Pd/P)^((gamma-1)/gamma)))^(1/2);
Pss=P*(1-Xi^2*(1-(Pd/P)^((gamma-1)/gamma)))^(gamma/(gamma-1));
Tss=T*(Pss/P)^((gamma-1)/gamma);
Td=T*(Pd/P)^((gamma-1)/gamma);
end

function [Pss , Tss, Uss, Pu , Tu]=Outlet_simplemodel(P,T,Xi, Amean,gamma,R,m)
Tu=T;
fun2= @(Pu) Pu/(R*Tu)*(1-Xi^2*(1-(P/Pu)^((gamma-1)/gamma)))^(1/(gamma-1))*Xi*Amean*(((2*gamma*R*Tu)/(gamma-1))*(1-(P/Pu)^((gamma-1)/gamma)))^(1/2)-m;
Purange=[P 6];
Pu=fzero(fun2,Purange);
Uss=((2*gamma*R*Tu)/(gamma-1)*Xi^2*(1-(P/Pu)^((gamma-1)/gamma)))^(1/2);
Pss=Pu*(1-Xi^2*(1-(P/Pu)^((gamma-1)/gamma)))^(gamma/(gamma-1));
Tss=Tu*(Pss/Pu)^((gamma-1)/gamma);
end
end 
