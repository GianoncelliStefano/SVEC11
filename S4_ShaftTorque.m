function [R_x,R_xg,R_y,R_yg,R_rot,C_rot,C_bronz,C_shaft,C_press,C_shaft_mean] = S4_ShaftTorque(C1_pal_rot,C1_pal_P,R1_x,R1_y,n_pal,N_vano,d_hub,f_b,Fg,zeta,rot_dir,c)
% Torque and forces acting on shaft due to a single vane are summed up to obtain the overall quantities
% INPUT
% C1_pal_rot [Nm]  : array of torque acting on the rotor due to a single vane
% C1_pal_P [Nm]    : array of torque acting on a single vane due to indicated pressure (necessary to indicated power calculation)
% R1_x [N]         : x-component of shaft reaction force due to a single vane
% R1_y [N]         : y-component of shaft reaction force due to a single vane
% n_pal [-]        : number of vanes
% N_vano [-]       : grid discretization points
% d_hub [m]        : shaft diameter
% f_b  [-]         : bushing friction coefficient
% Fg [N]           : global weight force  
% zeta [°]         : angles of x axis with respect to X
% rot_dir [-]      : direction of rotation     1: clockwise       2: counterclockwise
% c  [-]           : stator geometry     1: circular   2: elliptical stator;
%
% OUTPUT
% R_x  [N]         : overall x-component of shaft reaction force
% R_xg [N]         : x-component of the total resultant on the rotor accounting for the rotor-shaft weight force
% R_y  [N]         : overall y-component of shaft reaction force
% R_yg  [N]        : y-component of the total resultant on the rotor accounting for the rotor-shaft weight force
% R_rot [N]        : magnitude of shaft reaction force 
% C_rot  [Nm]      : overall torque acting on shaft due to reaction forces 
% C_bronz  [Nm]    : absolute value of torque due to bushing friction
% C_shaft  [Nm]    : overall mechanical torque acting on shaft, as function of vane angular position
% C_press  [Nm]    : overall torque acting on shaft due to pressure forces as function of vane angular position
%
% HISTORY: 
% Franzetti-Persico (june 2019) have removed the usage of s_bush variable that was used
% in Vallone to account for correct sign of C_bronz in case of
% compressor/expander.
    %% COMPUTATIONS %%
    % sum of effects on all vanes
    [C_rot]   = S4_AllVanes(C1_pal_rot,n_pal,N_vano);  % overall torque acting on shaft due to reaction forces [Nm]
    [C_press] = S4_AllVanes(C1_pal_P,n_pal,N_vano);    % overall torque acting on shaft due to pressure forces [Nm]
    [R_x]     = S4_AllVanes(R1_x,n_pal,N_vano);        % overall x-component of shaft reaction force [N]
    [R_y]     = S4_AllVanes(R1_y,n_pal,N_vano);        % overall y-component of shaft reaction force [N]

    if c == 1
        R_xg = R_x - Fg*sin(zeta);           % x-component of the total resultant on the rotor accounting for the rotor-shaft weight force
        R_yg = R_y - rot_dir*Fg*cos(zeta);   % y-component of the total resultant on the rotor accounting for the rotor-shaft weight force
    else
        % TO DO for c == 2 case;
        R_xg = R_x;
        R_yg = R_y;
    end
    
    % Total torque
    R_rot    = sqrt(R_xg.^2+R_yg.^2);      % magnitude of shaft reaction force [N] 
    C_bronz  = d_hub/2*R_rot*f_b;         % absolute value of torque due to bushing friction [Nm]
    C_shaft  = C_rot + C_bronz;   % total shaft torque [Nm]
    C_shaft_mean=mean(C_shaft); %PROVA_RICCARDO
end

%% Subfunction S4_AllVane
% this function takes as input a quantity referred to a single vane, and
% computes the corresponding overall quantity
function [C_tot] = S4_AllVanes(C,n_pal,N_vano)
    C     = reshape(C,N_vano,n_pal);
    C1    = sum(C,2);
    C1    = C1';
    C_tot = repmat(C1,1,n_pal);
end

