function [VB,VI] = S6_ApplicationPointForces(Y_F,HI, AV_i, d, th_tilt)
% This function computes the 1D geometry segments used in stress analysis.
% INPUT
% Y_F [m]        : y-coord of center of curvature of tip
% HI [m]         : tip excursion
% AV_i [m]       : from A to bottom of vane segment
% d  [m]         : rotor diameter
% th_tilt [rad]  : vane tilt angle
%
% OUTPUT
% VB: [m]  : segment from vane bottom to vane-rotor interseption (point B)
% VI: [m]  :  vane axis length

    %% 2D VANE GEOMETRY CALCULATIONS %%
    r = 0.5*d;                  % rotor radius [m]
    VB = r*cos(th_tilt) - AV_i; % VB segment [m]
    VI = Y_F + HI;              % VI segment [m]
    clear r

end
