function S7_PlotMechanics(F_m,F_v,F_t,T_m,T_v,T_t,F_cor,F_centr_i,F_centr_j,F_iz_i,conf,th_conf,Npt)
% This function create the following plots, as function of the angle on the vane (0:2pi)
% -reaction forces
% -friction forces
% -vane configuration
%
% INPUT: 
% F_m [N]          : reaction force in rotor-vane contact point
% F_v [N]          : reaction force at vane bottom
% F_t [N]          : reaction force at vane tip
% T_m [N]          : friction force in rotor-vane contact point (generated by F_mA) 
% T_v [N]          : friction force at vane bottom (generated by F_vA) 
% T_t [N]          : friction force at vane tip (generated by F_tA) 
% F_cor [N]        : Corolis force (perpendicular to vane)
% F_centr_i [N]    : i-comp. centrifugal force
% F_centr_j [N]    : j-comp. centrifugal force
% F_iz_i [N]       : inertial force (along the vane)
% conf [-]         : array of vane configurations
% Npt      [-]     : number of grid discretization points in [0:2pi]

    %% DEFINITIONS %%
    x = linspace(0,360,Npt);
    
    %% CONVERSION %%
    th_conf = rad2deg(th_conf);    % [red -> deg]
    
    %% PLOT REACTION FORCES (F_m F_v F_t) %%
    figure('Name','Fig_5 Vane forces');
    subplot (3,1,1);
    plot(x,F_m,'b',x,F_v,'r',x,F_t,'g'); hold on;
    SX_vline(th_conf,'k:');
    %options
    grid on, box on;
    xlim([0 360]);
    title('Reaction Forces');
    xlabel('Theta [deg]');
    ylabel('Force [N]');
    legend('F_m','F_v','F_t');

    %% PLOT FRICTION FORCES (T_m; T_v; T_t) %% 
    subplot (3,1,2);
    plot(x,T_m,'b',x,T_v,'r',x,T_t,'g'), hold on;
    SX_vline(th_conf,'k:');
    %options
    grid on, box on;
    xlim([0 360]);
    title('Friction Forces');
    xlabel('Theta [deg]');
    ylabel('Force [N]');
    legend('T_m','T_v','T_t');

    %% PLOT CONFIGURATIONS %%
    subplot (3,1,3);     
    plot(x,conf,'o'), hold on;
    SX_vline(th_conf,'k:');
    %options
    grid on, box on; 
    xlim([0 360]);
    title('Vane Configurations');
    xlabel('Theta [deg]');
    ylabel('Conf[-]');

    %% PLOT TANGENTIAL FORCES (F_cor; F_centr_j) %%
    figure('Name','Fig_6 Vane body forces');
    subplot (2,1,1);
    plot(x,F_cor,'b',x,repmat(F_centr_j,1,size(x,2)),'r');
    %options
    grid on, box on
    xlim([0 360]);
    ylim('auto');
    title('Tangential Forces');
    xlabel('Theta [deg]');
    ylabel('Force [N]');
    legend('Fcor','Fcentr_j');

    %% PLOT RADIAL FORCES (F_centr_i; F_iz_i) %%
    subplot (2,1,2);
    plot(x,F_centr_i,'b',x,F_iz_i,'r');
    %options
    grid on, box on;
    xlim([0 360]);
    title('Radial Forces');
    xlabel('Theta [deg]');
    ylabel('Force [N]');
    legend('Fcentr_i','Fiz_i');

end
