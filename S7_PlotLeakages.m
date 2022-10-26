function S7_PlotLeakages(theta_vane,theta_SucOpen,theta_SucClose,theta_DisOpen,theta_DisClose,Flow_RS,Flow_VS,Flow_PE,Gamma)
% This function creates the plot of the three different leakages
% implemented in SVEC
% INPUT:
% theta [rad]      : array of discretized angular position of trailing vane [-gamma:2pi]


% Gamma [rad]      : angular extention of a cell (cavity between 2 vanes)
% DPLT [bool]      : display plot       0: Disable        1: Enable
    %% DEFINITION %%
   
    %% CONVERSION %%
    %Gamma       = rad2deg(Gamma);      % [rad -> deg]
    %theta       = rad2deg(theta);      % [rad -> deg]
    %x_deg       = rad2deg(x_rad);      % [rad -> deg]
    
    %% PLOT LEAKAGES %%

    %Yrange
    YY_plus= max(Flow_RS)*1.5;
    YY_minus=-max(Flow_RS)*1.5;
         figure('Name','Fig_8 Leakages');
         plot(theta_vane,Flow_RS,'g');
         hold on
         grid on
         plot(theta_vane,Flow_PE,'r');
         plot(theta_vane,Flow_VS,'b');
         plot([theta_SucClose theta_SucClose],[YY_minus YY_plus],'k--')
         plot([theta_DisOpen-Gamma theta_DisOpen-Gamma],[YY_minus YY_plus],'k--')
         plot([theta_DisOpen theta_DisOpen],[YY_minus YY_plus],'k--')
         plot([theta_DisClose theta_DisClose],[YY_minus YY_plus],'k--')
         title('Axial leakages');
         text(theta_SucClose-0.4, YY_minus,'SUC_C','Color','black','FontSize',16);
         text(theta_DisOpen-Gamma-1, YY_plus,'DIS_O-Gamma','Color','black','FontSize',16);
         text(theta_DisOpen-0.5, YY_minus,'DIS_O','Color','black','FontSize',16);
         text(theta_DisClose, YY_minus,'DIS_C','Color','black','FontSize',16);
         xlabel('Trailing vane angle [deg]');
         ylabel('leakage [kg/s]');
         legend('rotor-stator','rotor-end-palte','vane','Location','northwest');
         hold off
    
    
    
end


