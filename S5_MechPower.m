function [P_mecc,P_ind,C_shaft] = S5_MechPower(C_shaft,C_press,omega,Pind_PV,toll_d)
% This function computes the istantaneous mechanical and indicated power
% as function of vane angular position.
% Then, the mean mechanical and indicated power are calculated
% INPUT
% C_shaft [Nm]   : overall mechanical torque acting on shaft, as function of vane angular position
% C_press [Nm]   : overall torque acting on shaft due to pressure forces as function of vane angular position
% omega [rad/s]  : vane angular velocity
% P_ind_pv [W]   : indicated power computed buy integration of PV diagram
% toll_d [adim]  : numerical tolerance
%
% OUTPUT 
% P_mecc [W]     : mean shaft mechanical power 
% P_ind [W]      : mean shaft indicated power - power used for compression only

    %% POWER COMPUTATION %%
    % istantaneous power (function of vane angular position)
    P_mecc_ist = C_shaft*omega;  % istantaneous shaft mechanical power as function of vane angular position [W]
    P_ind_ist  = C_press*omega;  % istantaneous shaft indicated power as function of vane angular position [W]
    
    % indicated power
    P_mecc = mean(P_mecc_ist);   % mean shaft mechanical power 
    P_ind  = mean(P_ind_ist);    % mean shaft indicated power [W] 

    %% CHECKS %%
    errP  = abs((P_ind-Pind_PV)/Pind_PV);
    if errP > toll_d
        warning('S4_Power:IndicatedPower','Indicated power do not match. Estimated relative error: %2.4f %%', errP*100)
        SX_Logfile ('w',{strrep(lastwarn,'%','%%')});
    end, clear errP  
end