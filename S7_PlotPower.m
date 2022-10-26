function S7_PlotPower (C1_pal_rot,C1_pal_P,C1_tip,C1_iz,C_izG,Npt)
% This function create the plot of different torque contribution, as function of the angle on the vane (0:2pi)
%
% INPUT:
% C1_pal_rot [Nm]  : array of torque acting on the rotor due to a single vane
% C1_pal_P [Nm]    : array of torque acting on a single vane due to indicated pressure (necessary to indicated power calculation)
% C1_tip [Nm]      : array of torque acting on a single vane due to tip friction
% C1_iz [Nm]       : array of torque acting on a single vane due to inertial effects
% C_izG [Nm]       : inertial torque
% Npt [-]          : number of grid discretization points in [0:2pi]

    %% DEFINITIONS %%
    x1 = linspace(0,360,Npt);
    x2 = x1(1:end-1);
    
    %% PLOT INERTIAL MOMENT %%
    % figure, hold on,grid on, box on;
    % plot(x1,C_izG);
    % xlim([0 360]);
    % title('Inertia Moment');
    % xlabel('theta [deg]');
    % ylabel('C [N m]');

    %% PLOT TORQUE %%
    figure('Name','Fig_7 Mech torque');
    p1 = plot(x2,C1_pal_rot,x2,C1_pal_P,x2,C1_iz,x2,-C1_tip);
    %options
    grid on, box on;
    xlim([0 360]);
    title('Mechanical Torque');
    xlabel('Theta [deg]');
    ylabel('Torque [Nm]');
    legend('T_{Tot Shaft}','T_{Fpress}','T_{Finert}','T_{Ftip}','AutoUpdate','off','Location','best');
    uistack(p1(1),'top');

end