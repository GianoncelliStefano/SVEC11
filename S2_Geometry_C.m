function [Rth,PSI,BI,HI,UPSILON,UPSILON_tip_1,UPSILON_tip_2,ME,NL,EP_arch,PL_arch,BM,BN,AX_i,AY_i,AP_i,AP_j,EI,IL,AI,AI_1D,AR_i,AR_j,AS_i,AS_j,AV_i,AM_i,R_v,VN_i,VM_i,sgn]...
    = S2_Geometry_C(D,d,theta_v,theta_vane,th_tilt,Sigma,r_tip,XI,b,s,DELTA,LAMBDA,toll_t,Y_F)
%this function solves the contact point between vane and stator through an
%analytical way. Then real excursion and real tip reaction inclination are
%evaluated. Finally the surfaces exposed to working fluid and forces arms
%are calculated
%
%INPUT
% D [m]                  : stator diameter
% d [m]                  : rotor diameter 
% theta_vane [rad]       : vane angular positions between (0:2pi)
% th_tilt [rad]          : vane tilt angle
% r_tip [m]              : tip radius
% b [m]                  : center of tip curvature offset
% s [m]                  : vane thickness
% Sigma, DELTA [rad]     : characteristic angles of vane
% XI,LAMBDA theta_v        
% toll_t [-]             : teorical tollerance
% Y_F [m]                : y coordinates of center of curvature of tip 
%
%OUTPUT
% Rth [m]                : stator radius respect to rotor center O
% PSI [rad]              : angle that describes the position of the contact point along the vane tip
% BI [m]                 : real vane excursion measured on vane axis
% HI [m]                 : vane tip excursion                            
% UPSILON [rad]          : tip reaction force angle
% UPSILON_tip_1 [rad]    : tip pressure force angle - suction side
% UPSILON_tip_2 [rad]    : tip pressure force angle - pressure side
% ME, NL, [m]            : gas-exposed vane surfaces (lateral surfaces)
% EP_arch, PL_arch [m]   : gas-exposed vane surfaces (on vane tip)
% AX_i [m]               : arm of exposed pressure force (1)
% AY_i [m]               : arm of exposed pressure force (2)
% AP_j,AP_i [m]          : arms of tip reaction forces
% AR_j, AR_i [m]         : arms of tip pressure force (1)
% AS_j, AS_i [m]         : arms of tip pressure force (2)
% AV_i [m]               : rotor bottom reaction force (Fv) arm
% AM_i [m]               : rotor top reaction force (Fm) arm
% BM, BN [m]             : rotor-vane chord
% EI, IL [m]             : tip-vane chord
% AI [m]                 : geometrical quantity for 2d vane
% AI_1D [m]              : geometrical quantity for 1d vane
% VN_i [m]               : vane side surfaces inside the rotor slot - suction side
% VM_i [m]               : vane side surfaces inside the rotor slot - pressure side
% sgn [-]                : vane tilt direction
%
%NOTES
% see D Zanchi thesis for references and pictures (Chapter 2)
% (1) refers to LOW pressure side, (2) refers to HIGH pressure side
% [compression mode]
% HISTORY: V10.1 DeFranco_Genoni_Gianoncelli: AV_i modified to adjust areas evaluation in some peculiar cases

    %% DECLARATIONS %%
    % tilted vane discriminant factor
    if th_tilt >= 0    
        sgn =  1;
    else
        sgn = -1;
    end
    
    r          = d/2;                                                    % rotor radius  [m]
    R          = D/2;                                                    % stator radius [m]
    e          = R-r;                                                    % eccentricity  [m]
    Rth        = sqrt(R^2 - (e*sin(theta_vane)).^2) - e*cos(theta_vane); % stator radius respect to rotor center O [m]
    r_v        = abs(r*sin(th_tilt));                                    % corresponds to OA [mm]
    eta        = sgn*(pi/2- abs(th_tilt));                               % AOB angle [rad]
    alpha      = theta_vane - eta;
    q          = theta_vane + th_tilt;
    
    %% SURFACE CALCULATIONS %%
    % calculation of contact point between blade tip and stator
    R_v = -(b+r*sin(theta_vane).*cos(q)-e*sin(q)-r*cos(theta_vane).*sin(q));
    PSI = asin(R_v/(R-r_tip))- theta_v;
    
    % real excursion calculation [m]
    HI    = r_tip*sin(XI);                                    
    AH    = sqrt((R - r_tip)^2 - R_v.^2) +e*sin(alpha)*sgn;
    AI    = AH + HI;
    AB    = r*cos(th_tilt);
    BI    = AI - AB;                                                
    AI_1D = sgn*e*sin(alpha) + sqrt(R^2 - (r_v + e*cos(alpha)).^2); % former AC segment used for 1D vane
    
    % geometrical quantities useful for exposed surfaces computation [m]
    BH   = BI - HI;
    BM_i = (r*cos(asin((r_v - sgn*0.5*s)/r)) - AB);
    BN_i = (-r*cos(asin((r_v + sgn*0.5*s)/r)) + AB);
    BM   = sqrt(BM_i^2 + (0.5*s)^2);
    BN   = sqrt(BN_i^2 + (0.5*s)^2);    
    ED   = r_tip*cos(DELTA);
    LQ   = r_tip*sin(LAMBDA);
   
    % exposed surfaces [m]
    ME      = BH + ED - BM_i;               % vane pressure force surface (2)
    NL      = BH + LQ + BN_i;               % vane pressure force surface (1)
    EP_arch = (Sigma + PSI).*r_tip;         % tip pressure force surface (1)
    PL_arch = (Sigma - PSI).*r_tip;         % tip pressure force surface (2)    
    EI      = sqrt((HI-ED)^2 + (0.5*s)^2);  % chord of arc EI
    IL      = sqrt((HI-LQ)^2 + (0.5*s)^2);  % chord of arh IL
    clear r R e r_v eta alpha q ED LQ
    
    %% ARMS & FORSE INCLINATION CALCULATIONS %%
    % vane pressure forces arm [m]
    AX_i = AB + 0.5*NL - BN_i; % vane pressure force arm (1)
    AY_i = AB + 0.5*ME + BM_i; % vane pressure force arm (2)
    
    UPSILON = theta_v + PSI;           % tip reaction force inclination
    
    % tip pressure force arm [m] and inclination [rad]
    if s == 0
        UPSILON_tip_1 = UPSILON;                        %1D geometry
        UPSILON_tip_2 = UPSILON;                        %1D geometry
    else
        UPSILON_tip_1 = (UPSILON + 0.5*pi - LAMBDA)/2;  % tip pressure force inclination (1) %UPSILON_tip_1 = theta_v + Sigma - 0.5*(Sigma - PSI);
        UPSILON_tip_2 = (UPSILON - DELTA)/2;            % tip pressure force inclination (2) %UPSILON_tip_2 = theta_v - Sigma + 0.5*(Sigma + PSI);
    end

    AR_i          = AH + r_tip*cos(UPSILON_tip_1);      % tip pressure force arm (1)
    AR_j          = -b + r_tip*sin(UPSILON_tip_1);      % tip pressure force arm (1)
    AS_i          = AH + r_tip*cos(UPSILON_tip_2);      % tip pressure force arm (2)
    AS_j          = -b + r_tip*sin(UPSILON_tip_2);      % tip pressure force arm (2)
    
    % tip reaction force arm [m] and inclination [rad]
    AP_i    = AH + r_tip*cos(UPSILON); % tip reaction force arm
    AP_j    = r_tip*sin(UPSILON) - b;  % tip reaction force arm
    
    % rotor reaction forces arms [m]
    AM_i      = NaN(2, length(theta_vane));
    AM_i(1,:) = AB - BN_i + NL.*(NL<0);   % rotor top reaction force arm (1)
    AM_i(2,:) = AB + BM_i + ME.*(ME<0);   % rotor top reaction force arm (2)
    AM_i      = [AM_i; flipud(AM_i)];
    AV_i      = AH - Y_F;                  % rotor bottom reaction force arm
    
    % vane lateral surfaces inside the rotor [m]
    VN_i = AM_i(1,:) - AV_i;
    VM_i = AM_i(2,:) - AV_i;      
    clear AH AB BH BM_i BN_i theta_vane
    
    %% CHECKS %%
    % contact point check
    check_1 = abs(PSI) > Sigma;
    if sum(check_1) > 0
        warning('S2_Geometry_C:Slippage','cuspidale point slippage occurs! please modify r_tip and b');
        SX_Logfile ('w',{lastwarn});
    end
    
    % control on real excursion: it must lower than 1D excursion
    delta_AI = (AI_1D-AI);
    delta_AI(abs(delta_AI)<toll_t) = 0;
    check_2 = delta_AI < 0;
    if sum(check_2) > 0
        warning('S2_Geometry_C:contactpoint','stator-vane contact point falls outside the stator for some angular positions');
        SX_Logfile ('w',{lastwarn});
    end
    
    % control on tip pressure force inclination
    UPSILON_tip_1_check = theta_v + Sigma - 0.5*(Sigma - PSI);
    UPSILON_tip_2_check = theta_v - Sigma + 0.5*(Sigma + PSI);
    check_3 = abs(UPSILON_tip_1 - UPSILON_tip_1_check) > toll_t;
    check_4 = abs(UPSILON_tip_2 - UPSILON_tip_2_check) > toll_t;
    if sum(check_3)>0 || sum(check_4)>0
        warning('S2_Geometry_C:wrong tip pressure force inclinations');
        SX_Logfile ('w',{lastwarn});
    end

    clear check_1 check_2 check_3 check_4 delta_AI UPSILON_tip_1_check UPSILON_tip_2_check
end