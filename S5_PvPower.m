function [Wtot,Wis_g,Wis_gl,Wis_T,Pind_PV] = S5_PvPower(W_cell,p_suc,p_out,V_suc,V_dis,rpm,n_van,c,R_g,c_v,c_l,T_suc,m_g,m_inj_nzl,V_inj_nzl)
% This function computes the work and power derived from thermodynamic process
% INPUT
% W_cell [J]     : work done by the cell as close chamber
% p_suc [Pa]     : suction pressure
% p_out [Pa]     : end compression pressure
% V_suc [m3]     : cell volume vector during suction
% V_asp [m3]     : cell volume vector during discharge
% rpm [rpm]      : rotational speed
% n_van [-]      : number of vanes
% c [-]          : flag for gemetry selection (concentrical/elliptical)
% R_g [J/kgK]    : gas specific constant
% c_v [J/kgK]    : specific heat at constant volume of gas
% c_l [J/kgK]    : specific heat of liquid
% T_suc [K]      : initial compression temperature
% m_g [kg]       : mass of gas in a closed cell
% eps_m [-]      : oil/gas mass ratio
% rho_o [kg/m3]  : oil density
% m_inj_nzl [kg] : mass of oil injected by each nozzle
% V_inj_nzl [m3] : volume of oil injected by each nozzle
%
% OUTPUT
% Wtot [J]    : work done by the cell during a complete rotation
% Wis_g [J]   : gas-only isentropic closed-cell work
% Wis_gl [J]  : gas-liquid isentropic closed-cell work
% Wis_T [J]   : isothermal closed-cell work
% Pind_PV [W] : indicated power
%
% NOTE: - Closed cell means during the phase between suction closing and delivery opening
%       - With the signs implemented the two impulse works are just summed 
%          because discharge is positive and suction is negative.

    %% DEFINITIONS %%
    c_p     = c_v+R_g;                    % gas specific heat at costant pressure [J/kgK]
    gamma   = c_p/c_v;                    % specific heat ratio
    THETA   = (gamma-1)/gamma;            % adiabatic exponent [adim]
    beta    = p_out/p_suc;                % compression ratio [adim]
    epsc    = c_l/c_p;                    % liquid-gas specific heat ratio [adim]
    m_l_tot = sum(m_inj_nzl);             % total mass of injected oilo [kg]
    V_l_tot = sum(V_inj_nzl);             % total volume of injected oilo [m3]
    cgl     = c_p*(1+m_l_tot/mean(m_g)*epsc);   % equivalent spacific heat of gas-liquid mixture [J/kgK] 
    THETAgl = R_g/cgl;                    % equivalent exponent of adiabatic process of gas-liquid mixture [adim]
    deltap  = p_out-p_suc;                % pressure delta between suction and delivery [pa]
    clear W_suc Wdis
    
    %% WORK %%
    % pulsion work calculation
    W_suc = -p_suc*diff(V_suc);  % impulse work at suction   [J] 
    W_dis = -p_out*diff(V_dis);  % impulse work at discharge [J]
    DPV = sum(W_suc)+sum(W_dis); % total impulse work [J] 
    
    % others thermodynamic works
    Wtot   = W_cell+DPV;                                              % open system real work for a single cell [J]
    Wis_g  = mean(m_g)*c_p*T_suc*(beta^THETA-1);                            % gas-only isentropic compression work [J] 
    Wis_gl = (mean(m_g)*cgl*T_suc*(beta^THETAgl-1) + sum(V_l_tot*deltap));  % gas-liquid isentropic compression work [J]
    Wis_T  = mean(m_g)*T_suc*R_g*log(beta);                                 % isothermal compression work [J]
    clear  cp THETA beta epsc m_l_tot cgl THETAgl deltap
    
    %% INDICATED POWER COMPUTATION - thermodynamic process %%
    t_cycle = 60/(c*rpm);          % time employed for a complete turning [s]
    Pind_PV = Wtot*n_van/t_cycle;  % power elaborated by the machine during a complete rotation
    clear t_ciclo Area omega

end