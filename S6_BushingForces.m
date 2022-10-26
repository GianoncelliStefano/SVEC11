function [FB,kB,pB,vB,pvB] = S6_BushingForces(bB,R_rot,R_xg,R_yg,rpm,d_hub,zeta,rot_dir)
% Bushing mechanical performance calculation
%
% INPUT
% bB      [m]            : bushing length 
% R_rot   [N]            : magnitude of shaft reaction force
% R_xg    [N]            : x-component of the total resultant on the rotor (including weight force)
% R_yg    [N]            : y-component of the total resultant on the rotor (including weight force)
% rpm     [adim]         : rotational speed
% d_hub   [m]            : bushing diameter 
% zeta    [°]            : angles of x axis with respect to X
% rot_dir [-]            : direction of rotation     1: clockwise       2: counterclockwise

%
% OUTPUT
% FB  [N]                : maximum bushing load
% kB  [deg]              : bushing load direction
% pB  [N/mm^2]           : bushing specific load
% vB  [m/s]              : sliding velocity on bushings
% pvB [N/mm^2 • m/s]     : bushing Pressure-Volume Value

%% BUSHING MECHANICAL PERFORMANCES CALCULATIONS %%
FB    = max(R_rot)/2;                                    % Maximum bushing load [N]
kB    = (zeta + rot_dir*atan(R_yg./R_xg))*180/pi;        % Bushing load direction [deg]
pB    = FB/(bB*1000*d_hub*1000);                         % Bushing specific load [N/mm^2]
vB    = 2*pi*rpm/60*d_hub/2;                             % Sliding velocity on bushings [m/s]
pvB   = pB*vB;                                           % Bushing Pressure-Volume Value [N/mm^2 • m/s]

end