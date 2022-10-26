function S7_PlotGeometry(PSI,Sigma,chck_tip,Npt)
%this function plots the tip contact point variation during an entire cycle
% INPUT:
% PSI [rad]    : angle that describes the position of the contact point along the vane tip
% Sigma [rad]  : LFE/2 angle (half of the tip sweep angle)
% N_pt [-]     : number of corrected grid points [0:2pi]
% chck_tip     : boolean array - 1: contact between vane and stator    0: NO contact between vane and stator

    %% DEFINITION %%
    theta_deg = linspace(0,360,Npt);
    
    %% CONVERSION %%
    Sigma = rad2deg(Sigma);  % [rad _> deg]
    PSI   = rad2deg(PSI);    % [rad _> deg]
    
    %% PLOT TIP CONTACT POINT %%
    figure('Name','Fig_1 Vane contact');
    subplot (2,1,1);
    plot(theta_deg,PSI,'b'), hold on;
    SX_hline([Sigma,-Sigma],'r','Tip Contact Limit');
    %options
    grid on, box on;
    xlim([0 360]);
    title('Tip contact Point');
    xlabel('Theta [deg]');
    ylabel('Tip contact angle [deg]');
    legend('Tip Contact Point','Location','northeast');
    
    %% PLOT VANE-STATOR CONTACT POINT %%
    subplot (2,1,2);
    plot(theta_deg,chck_tip,'bo'), hold on;
  %  SX_hline(0,'r','no contact'), hold on;
  %  SX_hline(1,'r','contact'), hold on;
    %options
    grid on, box on;
    xlim([0 360]);
    ylim([-0.25 1.25]);
    title('Vane-stator contact check');
    xlabel('Theta [deg]');
    ylabel('Contact');
    set(gca,'YTick',[0,1]);
    set(gca,'yticklabel',{'NO'; 'YES'});
    %legend('Tip Contact Point','Location','northeast');
end

