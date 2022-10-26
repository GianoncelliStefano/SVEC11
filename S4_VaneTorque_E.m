function [C1_pal_rot,C1_pal_P,C1_tip,C1_iz] = S4_VaneTorque_E(Rth,d,x_g,l,w,beta,F_t,F_m,F_v,T_t,F_p,F_cor,toll_t)
% This function computes the different torque acting on the vane and the
% torque exchanged between vane and shaft
% INPUT
% Rth [m]    : stator "radius" (vector from rotor center to stator surface)
% d [m]      : rotor diameter
% x_g [m]    : center of mass position respect to vane bottom
% l [m]      : vane heigth
% w [m]      : vane excursion
% beta [rad] : tip reaction and friction force angle
% F_t [N]    : reaction force at vane tip
% F_m [N]    : reaction force in rotor-vane contact point
% F_v [N]    : reaction force at vane bottom
% T_t [N]    : tip friction
% F_p [N]    : pressure force
% F_cor [N]  : Coriolis force
% 
% OUTPUT:
% C1_pal_rot [Nm]      : array of torque acting on the rotor due to a single vane
% C1_pal_P [Nm]        : array of torque acting on a single vane due to indicated pressure (necessary to indicated power calculation)
% C1_tip [Nm]          : array of torque acting on a single vane due to tip friction
% C1_iz [Nm]           : array of torque acting on a single vane due to inertial effects
%
% NOTE
% the torque calculated is referred to the effect of a single vane

    %% ROTOR TORQUE %%
    %torque due to indicated effects
    C1_pal_rot = (F_m.*d/2 - F_v.*(Rth-l));             % Torque acting on the shaft due to a single vane (torque due to reaction forces) [Nm]
    
    %% VANE TORQUE %%
    C1_pal_P   = F_p.*(Rth-w/2);                        % Torque acting on a single vane due to the effect of indicated pressure [Nm]
    
    % Torque acting on a single vane due to non indicated effects:
    C1_tip     = (T_t.*cos(beta)-F_t.*sin(beta)).*Rth;  % Torque due to tip friction [Nm]
    C1_iz      = F_cor.*x_g;                            % Torque due to inertial effects (Coriolis force) [Nm]
   
    %% FINAL CALCULATIONS %%
    % variables calculated are transferred on an elliptical basis (0:pi) --> (0:2pi)
    C1_pal_rot(end) = [];
    C1_pal_P(end)   = [];
    C1_tip (end)    = [];   
    C1_iz (end)     = [];
    C1_pal_rot      = [C1_pal_rot C1_pal_rot];
    C1_pal_P        = [C1_pal_P C1_pal_P];
    C1_tip          = [C1_tip  C1_tip];   
    C1_iz           = [C1_iz C1_iz];
    
    % Total torque acting on a vane
    C1_pal_tot = C1_pal_P+C1_tip+C1_iz;   % Total torque acting on a vane [Nm] 
    
    %% CHECKS %%
    errC  = abs((C1_pal_rot - C1_pal_tot) ./ C1_pal_rot) ;
    chck = errC > toll_t;
    if sum(chck) >= 1
        warning('S4_VaneTorque_E:Torque', 'Torque acting on vane and on shaft do not match. Estimated relative error: %.3d %%', max(errC)*100)
        SX_Logfile('w',{strrep(lastwarn,'%','%%')});
    end,clear errC chck
    
    
end
