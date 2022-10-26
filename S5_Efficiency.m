function [eta_g,eta_gl,eta_gT,eta_vol,eta_mecc,L_sp_m,L_sp_i,E_sp_m] = S5_Efficiency(Wtot,Wis_g,Wis_gl,Wis_T,Pind_PV,P_mecc,m_gas_Id,m_gas,q_gas,m_inj_nzl,toll_d,process)
% Computation thermodynamic and mechanical efficiency
% INPUT
% Wtot [J]        : work done by the cell during a complete rotation
% Wtot [J]        : work done by the cell during a complete rotation
% Wis_g [J]       : gas-only isentropic closed-cell work
% Wis_gl [J]      : gas-liquid isentropic closed-cell work
% Wis_T [J]       : isothermal closed-cell work
% Pind_PV [W]     : indicated power% P_mecc [W]      : mean shaft mechanical power 
% m_gas_Id [kg/s] : ideal mass flow rate (no leakages)
% m_gas [kg/s]    : mass flow rate of gas
% q_gas [m3/s]    : vlumetric flow rate of gas
% m_inj_nzl [kg]  : mass of oil injected by each nozzle
% toll_d          : numerical tollerance
% process         : flag which select the process (compression/expansion)
%
% OUTPUT
% eta_g  [%]     : thermodynamic efficiency with respect to isentropic gas compression 
% eta_gl  [%]    : thermodynamic efficiency with respect to isentropic gas-liquid compression 
% eta_gT  [%]    : thermodynamic efficiency with respect to isothermal compression
% eta_vol [%]    : gas volumetric efficiency
% eta_mecc [%]   : mechanic efficiency
% L_sp_m [J/kg]  : specific mechanical massic work ( relative to mechanical power)
% L_sp_i [J/kg]  : specific indicated massic work ( relative to indicated power)
% E_sp_m [J/m3]  : specific mechanical volumetric work 
%
% NOTE: - Closed cell means during the phase between sunction closing and delivery opening

    %% EFFICIENCIES %%
    % compression efficiency
    if process == 1
    eta_g    = Wis_g/Wtot*100;       % efficincy respect to gas-only compression
    eta_gl   = Wis_gl/Wtot*100;      % efficiency respect to gas-liquid compression
    eta_gT   = Wis_T/Wtot*100;       % efficiency respect to isothermal compression
    eta_vol  = (m_gas)/m_gas_Id*100; % gas volumetric efficiency
    eta_mecc = Pind_PV/P_mecc*100;     % mechanical efficiency
    
    % expansion efficiencies
    elseif process == 2
       eta_g    = Wtot/Wis_g*100;
       eta_gl   = Wtot/Wis_gl*100;
       eta_gT   = Wtot/Wis_T*100;
       eta_vol  = m_gas_Id/m_gas*100;
       eta_mecc = P_mecc/Pind_PV*100;
    end

    %% SPECIFIC WORK and ENERGY [J/kg]
    L_sp_m = P_mecc/m_gas;   % specific mechanical massic work [J/kg]
    L_sp_i = Pind_PV/m_gas;    % specific indicated massic  work [J/kg]
    E_sp_m = P_mecc/q_gas;   % specific mechanical volumetric work [J/m3]
    
    %% CHECKS %%
    if eta_vol <0
        warning('S3_Efficiencies:Volumetric','Leakages effects are too high: all mass leaks. Chose more conservative leakage models')
        SX_Logfile ('w',{lastwarn})
    end
    if isempty(m_inj_nzl)
        errW = ((Wis_g-Wis_gl)/Wis_gl);
        if errW > toll_d
            warning('S3_Efficiencies:NoLiquidWork','No liquid present, but Gas-only compression work not corresponds to gas-liquid compression work.Estimated relative error: %2.4f %%',errW*100)
            SX_Logfile ('w',{lastwarn});
        end
    end
end