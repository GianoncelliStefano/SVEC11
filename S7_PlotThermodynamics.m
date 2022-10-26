function S7_PlotThermodynamics(theta,V_cell,V_comp,V_disI,p_cell,p_suc,p_out,T_g,T_mix,c_v,R_g,th_in,th_out,Gamma)
% This function creates the plots for the PV cycle and the temperature
% INPUT:
% theta [rad]      : array of discretized angular position of trailing vane [-gamma:2pi]
% V_cell [m^3]     : cell volumes vector during a complete vane revolution (-gamma:2pi)
% V_comp  [m^3]    : cell volumes vector during closed cel phase
% V_disI  [m^3]    : volume of first cell after delivery opening
% p_cell [Pa]      : pressure vector of the cell [-gamma:2pi]
% p_suc [Pa]       : inlet gas pressure
% p_out [Pa]       : outlet gas pressure
% T_g [K]          : gas temperature during closed chamber phase
% T_mix [K]        : adiabatic mixing temperature during closed chamber phase
% c_v [J/kgK]      : specific heat at constant volume of gas
% R_g [J/kgK]      : gas specific constant
% th_in [rad]      : suction close angle (takes into account if user wants to maximise suction volume)
% th_out [rad]     : delivery open angle (takes into account expansion ports flipping)
% Gamma [rad]      : angular extention of a cell (cavity between 2 vanes)
% DPLT [bool]      : display plot       0: Disable        1: Enable
    %% DEFINITION %%
    NptPlot                      = length(T_g);                                         % number of elements on X axis
    x_rad                        = linspace(th_in,th_out-Gamma,NptPlot);                % X axis in radiant
    [p_isoT,p_adb,T_adb,V_plot]  = idealgas(p_suc,p_out,V_comp,V_disI,T_g(1),c_v,R_g);  % thermodynamics of noble process
    
    %% CONVERSION %%
    Gamma       = rad2deg(Gamma);      % [rad -> deg]
    theta       = rad2deg(theta);      % [rad -> deg]
    x_deg       = rad2deg(x_rad);      % [rad -> deg]
    p_cell      = p_cell*1e-5;         % [Pa -> bar]
    p_isoT      = p_isoT*1e-5;         % [Pa -> bar]
    p_adb       = p_adb*1e-5;          % [Pa -> bar]
    V_cell      = V_cell*1e9;          % [m3 -> mm3]
    V_plot      = V_plot*1e9;          % [m3 -> mm3]
    T_g         = T_g-273.15;          % [K -> °C]
    T_mix       = T_mix-273.15;        % [K -> °C]
    T_adb       = T_adb-273.15;        % [K -> °C]
    
    %% PLOT PV CYCLE %%
    figure('Name','Fig_2 PV plot');
    p1 = plot(V_cell,p_cell,V_plot,p_adb,V_plot,p_isoT);
    %options
    grid on, box on;
    xlim([0 max(V_cell)]);
    title('Pressure-Volume cycle');
    xlabel('Cell Volume [mm^3]');
    ylabel('Cell Presure [bar]');
    legend({'PV_{gas+oil}','PV_{isoS, gas}','PV_{isoT, gas}'},'AutoUpdate','off','Location','northeast');
    uistack(p1(1),'top');
    
    %% PLOT VOLUME %%
    figure('Name','Fig_3 Thermodynamic ');
    subplot(3,1,1);
    plot(theta,V_cell); hold on;
    plot(0,0,'k:'); hold on;
    SX_vline([x_deg(1),x_deg(end)],'k:');
    %options
    grid on, box on;
    xlim([-Gamma, 360]);
    title('Volume');
    xlabel('Trailing vane angle [deg]');
    ylabel('Volume [mm3]');
    legend('Cell volume','Closed cell');
    
    %% PLOT PRESSURE %%
    subplot(3,1,2);
    p3 = plot(theta(1:end-1),p_cell(1:end-1),x_deg,p_adb(1:end-1),x_deg,p_isoT(1:end-1)); hold on;
    SX_vline([x_deg(1),x_deg(end)],'k:');
    %options
    grid on, box on;
    xlim([-Gamma, 360]);
    title('Cell pressure');
    xlabel('Trailing vane angle [deg]');
    ylabel('Pressure [bar]');
    legend('P_{gas+oil}','P_{isoS, gas}','P_{isoT, gas}','AutoUpdate','off');
    uistack(p3(1),'top');
    
   %% PLOT TEMPERATURE %%
    subplot(3,1,3);
    p4 = plot(x_deg,T_g,x_deg,T_adb,x_deg,T_mix); hold on;
    SX_vline([x_deg(1),x_deg(end)],'k:');
    %options
    grid on, box on;
    xlim([-Gamma, 360]);
    title('Cell temperature');
    xlabel('Trailing vane angle [deg]');
    ylabel('Temperature [°C]');
    legend('T_{gas}','T_{isoS, gas}','T_{adb mix, gas+oil}','AutoUpdate','off');
    uistack(p4(1),'top');
end


function [p_isoT,p_adb,T_adb,V_plot] = idealgas(p_suc,p_out,V_comp,V_disI,T_gIN,c_v,R_g)
% This function computes the pressure, temperature and volume 
% of a isothermal and adiabatic process
% INPUT
% p_suc [Pa]      : inlet gas pressure
% p_out [Pa]      : outlet gas pressure
% V_comp  [m^3]   : cell volumes vector during closed cel phase
% V_disI  [m^3]   : first cell volume of delivery process
% T_gIN [K]       : inlet gas temperature
% c_v [J/kgK]     : specific heat at constant volume of gas
% R_g [J/kgK]     : gas specific constant
%
% OUTPUT
% p_isoT          : pressure of an isotermal process for ideal gas
% p_adb [pa]      : pressure of an adiabatic reversible process for ideal gas
% T_adb [K]       : temperature of an adiabatic reversible process for ideal gas
% V_plot [pa]     : vector of closed chamber volume  and the volume of the first cell after delivery open

    k               = 1+R_g/c_v;                          % politropic exponent for a adiabatic/reversible transformation
    Nptcls          = length(V_comp);                     % closed chamber discretization points
    p_isoT          = NaN(1,Nptcls+1);                    % prellocation
    p_adb           = NaN(1,Nptcls+1,1);                  % prellocation
    T_adb           = NaN(1,Nptcls+1,1);                  % prellocation
    p_isoT(1:end-1) = p_suc*V_comp(1)./V_comp;            % pressure of an isotermal process for ideal gas
    p_adb(1:end-1)  = p_suc*(V_comp(1)./V_comp).^k;       % pressure of an adiabatic reversible process for ideal gas
    T_adb(1:end)    = (p_adb./p_suc).^(1-1/k).*(T_gIN);   % ideal gas temperature of an adiabatic process
    p_isoT(end)     = p_out;                              % to close the diagram
    p_adb(end)      = p_out;                              % to close the diagram
    T_adb(end)      = [];                                 % 
    V_plot          = [V_comp, V_disI];                   % cell volume vector useful for some plot
    
end


