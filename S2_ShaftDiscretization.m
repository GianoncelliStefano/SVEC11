function [xsv] = S2_ShaftDiscretization(bB,d1,L,shm)
% This function perform the discretization of the the shaft along the axis
% INPUT
% bB  [m]    : Bushing length
% d1  [m]    : Rotor end-bushing distance
% L   [m]    : Rotor length
% shm [-]    : Rotor+shaft axial discretization points
%
% OUTPUT
% xsv [m]    : Rotor+shaft axial length discretization

   %% CALCULATION %%
   Lsh = bB + d1 + L + d1 + bB;     %[m] Total length for the rotor-shaft assembly
   xsv = linspace(eps,Lsh,shm);
