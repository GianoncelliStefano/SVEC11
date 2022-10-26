function S7_PlotStress(posMaxMf,N_posMaxMf,T_posMaxMf,Mf_posMaxMf,sVB,sBI,stp,pos_contct)
% This function create the plots of the internal axial, shear and bending
% action for the most strained angular position (where the maximum banding
% moment has been found). Then perform a 3d plot of the above stress for an
% entire vane revolution
% INPUT
% posMaxMf [rad]       : discretized angular position where maxMf has been found
% N_posMaxMf [N]       : internal action N for the angular position and vane abscissa where maxMf has been found
% T_posMaxMf [T]       : internal action T for the angular position and vane abscissa where maxMf has been found
% Mf_posMaxMf [Nm]     : internal action Mf for the angular position and vane abscissa where maxMf has been found
% sVB [-]              : cell-array of DB discretized lenght 
% sBI [-]              : cell-array of BC discretized lenght
% stp [m]              : vane discretization step for stress analysis
% s_contact [m]        : vane discretization step for stress analysis
% Npt [-]              : number of grid discretization points in [0:2pi]

    %% DEFINITIONS %%
    x1 = 1:length(sVB{1,posMaxMf});
    x2 = length(sVB{1,posMaxMf}):length(sVB{1,posMaxMf})+length(sBI{1,posMaxMf})-1;
    x  = [x1 x2];                       % vane abscissa [adim]
        
    %% CONVERSION %%
    x          = x*stp*1e3;             % [adim -> m -> mm] vane abscissa
    pos_contct = pos_contct *stp*1e3;   % [adim -> m -> mm] vane abscissa of vane-rotor contact point
    
    %% PLOT BENDING ACTION Mf %%
    figure('Name','Fig_8 Vane internal actions');
    subplot(3,1,1);
    plot(x,Mf_posMaxMf,'o'), hold on;
    SX_vline(pos_contct(posMaxMf), 'k--');
    % options
    grid on, box on;
    title('Bending action - Mb');
    xlabel('Vane abscissa [mm]');
    ylabel('Bending moment [Nm]');
    xlim([0,x(end)]);
    
    %% PLOT AXIAL ACTION N %%
    subplot(3,1,2)
    plot(x,N_posMaxMf,'o'), hold on;
    SX_vline(pos_contct(posMaxMf), 'k--');
    % options
    grid on, box on;
    title('Axial action - N');
    xlabel('Vane abscissa [mm]');
    ylabel('Force [N]');
    xlim([0,x(end)]);

    %% PLOT SHEAR ACTION T %%
    subplot(3,1,3)
    plot(x,T_posMaxMf,'o'), hold on;
    SX_vline(pos_contct(posMaxMf), 'k--');
    %options
    grid on, box on;
    title('Shear action - T');
    xlabel('Vane abscissa [mm]');
    ylabel('Force [N]');
    xlim([0,x(end)]);

end