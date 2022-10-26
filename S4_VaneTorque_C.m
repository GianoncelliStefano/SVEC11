function [C1_pal_rot,C1_pal_P,C1_tip,C1_iz] = S4_VaneTorque_C(AG,b_gx,s,rv,F_t,T_t,F_m,T_m,F_v,T_v,F_cor,F_centr_i,F_centr_j,Fw_i,Fw_j,F_iz_i,Fp_1,Fp_2,Fp_1_tip,Fp_2_tip,Fp_cava,Fp_1_cava,Fp_2_cava,C_izG,UPSILON,UPSILON_tip_1,UPSILON_tip_2,AX_i,AY_i,AP_i,AP_j,AR_i,AR_j,AS_i,AS_j,AV_i,AM_i_new,Ay_i,Ax_i,sgn,toll_t,s1,s2)
% This function computes the different torque acting on the vane and the
% torque exchanged between vane and shaft
% INPUT
% AG [m]               : AG segment
% b_gx [m]             : center of mass offset respect to vane axis
% s [m]                : vane thickness
% rv [m]               : raggio circonferenza definita dal punto di tangenza A
% F_t [N]              : reaction force at vane tip
% T_t                  : friction force at vane tip (generated by F_t) 
% F_m [N]              : reaction force in rotor-vane contact point
% T_m                  : friction force in rotor-vane contact point (generated by F_m) 
% F_v [N]              : reaction force at vane bottom
% T_v                  : friction force at vane bottom (generated by F_vA) 
% F_p [N]              : pressure force in exposed side
% F_cor [N]            : Corolis force (perpendicular to vane)
% F_iz_i [N]           : inertial force (along the vane)
% F_centr_i [N]        : i-comp. centripetal force
% F_centr_j [N]        : j-comp. centripetal force
% Fw_i [N]             : i-comp. weigth force [N]
% Fw_j [N]             : j-comp. weigth force [N]
% Fp_1 [N]             : pressure force on exposed suction side
% Fp_2 [N]             : pressure force on exposed pressure side
% Fp_1_tip [N]         : tip pressure force on suction side
% Fp_2_tip [N]         : tip pressure force on pressure side
% Fp_cava [N]          : pressure force on bottom of the vane
% Fp_1_cava [N]        : pressure force on concealed suction side of the vane
% Fp_2_cava [N]        : pressure force on concealed pressure side of the vane
% C_izG [Nm]           : inertial torque
% UPSILON [rad]        : tip reaction force angle
% UPSILON_tip_1 [rad]  : tip pressure force angle - suction side
% UPSILON_tip_2 [rad]  : tip pressure force angle - pressure side
% AX_i [m]             : arm of pressure force - exposed pressure side (2)
% AY_i [m]             : arm of pressure force - exposed suction side (1)
% AP_i, AP_j [m]       : arms of tip reaction forces
% AR_i, AR_j [m]       : arms of tip pressure force - suction side (1)
% AS_i, AS_j [m]       : arms of tip pressure force - pressure side (2)    
% AV_i [m]             : arm of reaction force at vane bottom (F_v,T_v)
% AM_i_new [m]         : arm of F_m force, according to configuration 
% Ay_i [m]             : arms of Fp_1_cava 
% Ax_i [m]             : arms of Fp_2_cava
% sgn [-]              : vane tilt angle discriminant
% toll_t [adim]        : tolleranza teorica
% s1                   : sign of forces F_m
% s2                   : sign of forces  F_v
%
% OUTPUT
% C1_pal_rot [Nm]      : array of torque acting on the rotor due to a single vane
% C1_pal_P [Nm]        : array of torque acting on a single vane due to indicated pressure (necessary to indicated power calculation)
% C1_tip [Nm]          : array of torque acting on a single vane due to tip friction
% C1_iz [Nm]           : array of torque acting on a single vane due to inertial effects
%
% NOTE
% the torque calculated is referred to the effect of a single vane

    %% ROTOR TORQUE %%   
    % Torque acting on the shaft due to a single vane (torque due to reaction forces) [Nm]
    C1_pal_rot = F_m.*AM_i_new - F_v.*AV_i - sgn*T_m.*(rv + sgn*s1*0.5*s) - sgn*T_v.*(rv - sgn*s2*0.5*s) - (Fp_2_cava.*Ay_i - Fp_1_cava.*Ax_i - sgn*Fp_cava*rv); 
    
    %% VANE TORQUE %%
    % Torque acting on a single vane due to the effect of indicated pressure (torque due to pressure on exposed side & tip) [Nm]
    C1_pal_P   = -Fp_1.*AX_i + Fp_2.*AY_i + Fp_1_tip.*(-sin(UPSILON_tip_1).*AR_i + sgn*cos(UPSILON_tip_1).*(rv+sgn*AR_j)) +...  
                  Fp_2_tip.*(-sin(UPSILON_tip_2).*AS_i + cos(UPSILON_tip_2).*(sgn*rv+sgn*AS_j));
    
    % Torque acting on a single vane due to non indicated effects:
    % C1_tip  [Nm]   : torque due to tip friction
    % C1_iz   [Nm]   : torque due to inertial effects
    % C1_vane [Nm]   : torque due to pressure in the concealed side 
    C1_tip = F_t.*(-sin(UPSILON).*AP_i + sgn*cos(UPSILON).*(rv+sgn*AP_j)) + T_t.*(cos(UPSILON).*AP_i + sgn*sin(UPSILON).*(rv+sgn*AP_j));   
    C1_iz  = (sgn*(F_cor + F_centr_j) - Fw_j).*AG + (-F_centr_i + F_iz_i + Fw_i).*(rv-abs(b_gx)) + C_izG;                          
    %C1_vane =  Fp_2_cava.*Ay_i - Fp_1_cava.*Ax_i - sgn*Fp_cava*rv;
    
    %% FINAL CALCULATIONS %%
    % the last element is canceled since it is a repeated element equal to the first one
    C1_pal_rot(end)  = [];
    C1_pal_P(end)    = [];
    C1_tip (end)     = [];
    C1_iz(end)       = [];
    %C1_vane(end)     = [];
    
    % Total torque acting on a vane [Nm] 
    C1_pal_tot  = C1_pal_P  + C1_iz + C1_tip; % + C1_vane;
    
    %% CHECKS %%
    % Check that torque exchanged between vane and rotor are equal
    errC = abs((C1_pal_rot - C1_pal_tot) ./ C1_pal_rot);
    chck = errC > toll_t;
    if sum(chck) >= 1
        warning('S4_VaneTorque_C:Torque', 'Torque acting on vane and on shaft do not match. Estimated relative error : %2.4f %% ', max(errC)*100)
        SX_Logfile ('w',{strrep(lastwarn,'%','%%')});
    end, clear chck errC
    
end
