function [Npt,Npt_cell, Gamma, theta, pos_tang, theta_vane] = S2_Discretization(N_pt_i,n_van,c)
%Compute the number of the grid discretization points used in computation
%
% INPUT
% N_pt_i  [-]   : number of selected discretization point
% n_van   [-]   : number vanes
% c       [-]   : geometry flag
%   1               :circular
%   2               :elliptical
%
% OUTPUT
% Npt [-]            : number of grid discretization points in [0:2pi]
% Npt_cell [-]       : number of discretization points between 2 vanes
% Gamma [rad]        : angular extention of a cell (cavity between 2 vanes)
% theta [rad]        : array of discretized angular position of trailing vane [-gamma:2pi]
% theta_vane [rad]   : array of discretized angular position  [0:2pi]
% pos_tang [-]       : index of array theta for the discretised point closest to tangency  

    %% GRID DETERMINATION %%
    Npt_cell = ceil((N_pt_i-1)/n_van)*c; % number of points between 2 vanes in the grid
    Npt      = (Npt_cell*n_van)/c+1;     % number of corrected grid points [0:2pi]
    Gamma    = 2*pi/n_van;               % cavity angular extension [rad]

    % array of discretized angular position of trailing vane   
    theta = linspace(-Gamma,2*pi/c,Npt+Npt_cell); % (-gamma:2pi)  circular / % (-gamma: pi) elliptical
                                                  
    % tangency position in theta array  - because of discrretization, tangency may occour in points different from theta=0;
    [~,pos_tang] = min(abs(theta-0));  
    theta_vane   = theta(pos_tang:end);      % vane angular positions between (0:2pi)[rad]
    
    %% CHECKS %%
    cond_Npt = (Npt==N_pt_i);
    if cond_Npt==0
        warning('S2_Discretization:points','Number of discretization points changed. Discretization points currently used : %i ', Npt)
        SX_Logfile ('v',{lastwarn});
    end

end