function [F_tA,T_tA,F_mA,T_mA,F_vA,T_vA,F_cor,F_iz_i,F_centr_i,F_centr_j,Fw_i,Fw_j,Fp_1,Fp_2,Fp_1_tip,Fp_2_tip,Fp_cava,Fp_1_cava,Fp_2_cava,C_izG,AM_i_new,Ay_i,Ax_i,Mf_ax,conf,th_conf,s1,s2,chck_tip] =...
    S4_VaneDynamics_C(f_c,f_t,L,s,th_tilt,theta_vane,mass_vane,moment_inertia_vane,AG,b_gx,UPSILON,UPSILON_tip_1,UPSILON_tip_2,ME,NL,EP_arch,PL_arch,AX_i,AY_i,AP_j,AP_i,AR_j,AR_i,AS_j,AS_i,AV_i,AM_i,VN_i,VM_i,vAGi,aGcentr_i,aGcentr_j,aGi,aGcor_j,DomegaG,sgn,p_1,p_2,Npt,rot_dir,zeta,toll_t,fDBG)
% This function calculates the forces acting on the vane for each discretized 0-to-360 angular position
% for a cylindrical stator
% INPUT
% f_c [-]              : slot friction coefficient
% f_t [-]              : vane tip friction coefficient
% L [m]                : rotor length
% l [m]                : vane heigth
% s [m]                : vane thickness
% th_tilt [rad]        : vane tilt angle
% theta_vane [rad]     : vane angular position between [0 : 2pi]
% mass_vane [kg]       : vane mass
% moment_inertia_vane  : vane moment of inertia 
% AG [m]               : AG segment
% b_gx [m]             : center of mass offset respect to vane axis
% UPSILON [rad]        : tip reaction force angle
% UPSILON_tip_1 [rad]  : tip pressure force angle - suction side
% UPSILON_tip_2 [rad]  : tip pressure force angle - pressure side
% ME [m]               : pressure side exposed vane surface (2)
% NL [m]               : suction side exposed vane surface (1)
% EP_arch [m]          : pressure side tip arc length (2) 
% PL_arch [m]          : suction side tip arc length (1) 
% AX_i [m]             : arm of pressure force - exposed pressure side (2)
% AY_i [m]             : arm of pressure force - exposed suction side (1)
% AP_i, AP_j [m]       : arms of tip reaction forces
% AR_i, AR_j [m]       : arms of tip pressure force - suction side (1)
% AS_i, AS_j [m]       : arms of tip pressure force - pressure side (2)    
% AV_i [m]             : arm of reaction force at vane bottom (F_v,T_v)
% AM_i [m]             : arm of reaction force at vane bottom (F_m,T_m)
% VN_i [m]             : vane side surfaces inside the rotor slot - suction side
% VM_i [m]             : vane side surfaces inside the rotor slot - pressure side
% vAGi [m/s]           : i-comp. of center of mass velocity
% aGcentr_i [m/s2]     : i-comp. of center of mass centrifugal acceleration
% aGcentr_j [m/s2]     : j-comp. of center of mass centrifugal acceleration
% aGi_A     [m/s2]     : i-comp. of center. of mass inertial acceleration (analytical)
% aGcor_j   [m/s2]     : j-comp. of center of mass Coriolis acceleration
% DomegaG  [rad/s2]    : center of mass angular acceleration
% sgn [adim]           : vane tilt angle discriminant
% p_1 [pa]             : array of low pressure acting on vane in [0:2pi]
% p_2 [pa]             : array of high pressure acting on vane in [0:2pi]
% Npt [-]              : grid discretization points
% rot_dir              : direction of rotation     1: clockwise       2: counterclockwise
% zeta                 : angles of x axis with respect to X
% toll_t [-]           : theorical tolerance
% fDBG   [-]           : developer flag
%
% OUTPUT
% F_tA [N]       : reaction force at vane tip
% T_tA [N]       : friction force at vane tip (generated by F_tA) 
% F_mA [N]       : reaction force in rotor-vane contact point
% T_mA [N]       : friction force in rotor-vane contact point (generated by F_mA) 
% F_vA [N]       : reaction force at vane bottom
% T_vA [N]       : friction force at vane bottom (generated by F_vA) 
% F_p [N]        : pressure force in exposed side
% F_cor [N]      : Corolis force (perpendicular to vane)
% F_iz_i [N]     : inertial force (along the vane)
% F_centr_i [N]  : i-comp. centrifugal force
% F_centr_j [N]  : j-comp. centrifugal force
% Fw_i [N]       : i-comp. weigth force [N]
% Fw_j [N]       : j-comp. weigth force [N]
% Fp_1 [N]       : pressure force on exposed suction side
% Fp_2 [N]       : pressure force on exposed pressure side
% Fp_1_tip [N]   : tip pressure force on suction side
% Fp_2_tip [N]   : tip pressure force on pressure side
% Fp_cava [N]    : pressure force on bottom of the vane
% Fp_1_cava [N]  : pressure force on concealed suction side of the vane
% Fp_2_cava [N]  : pressure force on concealed pressure side of the vane
% C_izG [Nm]     : inertial torque
% AM_i_new [m]   : arm of F_m force, according to configuration 
% Ay_i [m]       : arms of Fp_1_cava 
% Ax_i [m]       : arms of Fp_2_cava
% Mf_ax [Nm]     : torque due to axial forces
% conf [-]       : array of vane configurations
% th_conf [rad]  : array of angle where there is a change in configuration
% s1 [-]         : sign of F_m
% s2 [-]         : sign of F_v
% chck_tip       : boolean array - 1: contact between vane and stator    0: NO contact between vane and stator

%
% NOTE:  1) i-compp. is along the vane
%        2) j-comp. is perpendicular to the vane
%        3) pressure and suction side are referred to compressor logic
% - Hp:  - vane are considered as omogeneous, without holes
%        - oil is not considered
%        - costant pressure distribution on suction side due to discarge holes
%        - linear pressure distribution on pressure side (p_1 to p_2 linear) 
% Gianoncelli-Genoni
%        4)In line 181, the last term of b2 differs from the
%        one presented in De Franco thesis. The right one is the present in the code,+ Fp_2_tip.*sin(UPSILON_tip_2). 
%        5)This two corrections might sound unuseful (s1*s1=1 ; s2*s2=1 always).
%        however, the implementation is finally coherent with references (see
%        material provided by Gianoncelli_Genoni)
% See A CORTESE thesis for more info

%% DECLARATION %%
    method = 'p_min';        % "De Franco 10.0" use 'p_min', 'p_med', 'p_max' to choose the pressure acting on the base of the vane
                             % For previous SVEC versions the pressure acting on the base of the vane is always 'p_min'

    %% PRELIMINARY CALCULATIONS %%
    Npt_v   = length (theta_vane);
    g_const = SX_Constant({'Gravity'});
    
    %% KNOWN FORCES CALCULATIONS %%
    % pressure forces [N]
    Fp_1         = NL.*p_1*L;      % pressure force on exposed suction side
    Fp_2         = ME.*p_2*L;      % pressure force on exposed pressure side
    Fp_1(Fp_1<0) = 0;
    Fp_2(Fp_2<0) = 0;
    Fp_1_tip  = PL_arch.*p_1*L;     % tip pressure force on suction side
    Fp_2_tip  = EP_arch.*p_2*L;     % tip pressure force on pressure side
   
    switch method
        case 'p_min'
            p_eq_slot = (p_1 + p_2)/2;      % mean pressure on concealed pressure side - see NOTE
            Fp_1_cava = VN_i.*p_1*L;        % pressure force on concealed suction side of the vane
            Fp_2_cava = VM_i.*p_eq_slot*L;  % pressure force on concealed pressure side of the vane
            Fp_cava   = s*p_1*L;            % pressure force on bottom of the vane
            
            % arms of Fp_1_cava and Fp_2_cava
            Ax_i = AV_i + 0.5*VN_i;                                % constant pressure distribution (1)
            Ay_i = AV_i + (p_1 + 2.*p_2)./(3.*(p_1 + p_2)).*VM_i;  % linear pressure distribution (2) - rectangular trapezium
        case 'p_med'
            p_med     = (p_1 + p_2)/2;
            p_eq_slot = (p_med + p_2)/2;      % mean pressure on concealed pressure side - see NOTE
            Fp_1_cava = VN_i.*p_1*L;        % pressure force on concealed suction side of the vane
            Fp_2_cava = VM_i.*p_eq_slot*L;  % pressure force on concealed pressure side of the vane
            Fp_cava   = s*p_med*L;            % pressure force on bottom of the vane
            
            % arms of Fp_1_cava and Fp_2_cava
            Ax_i = AV_i + 0.5*VN_i;                                % constant pressure distribution (1)
            Ay_i = AV_i + (p_med + 2.*p_2)./(3.*(p_med + p_2)).*VM_i;  % linear pressure distribution (2) - rectangular trapezium
        case 'p_max'
            Fp_1_cava = VN_i.*p_1*L;        % pressure force on concealed suction side of the vane
            Fp_2_cava = VM_i.*p_2*L;  % pressure force on concealed pressure side of the vane
            Fp_cava   = s*p_2*L;            % pressure force on bottom of the vane
            
            % arms of Fp_1_cava and Fp_2_cava
            Ax_i = AV_i + 0.5*VN_i;                                % constant pressure distribution (1)
            Ay_i = AV_i + 0.5*VM_i;  % linear pressure distribution (2) - rectangular trapezium
    end

    
    % Weight and inertia forces
    theta_kin = theta_vane+th_tilt;
    F_centr_i = mass_vane*aGcentr_i;              % i-comp. centrifugal force [N]
    F_centr_j = mass_vane*aGcentr_j;              % j-comp. centrifugal force [N]
    F_iz_i    = mass_vane*aGi;                    % inertial force (along the vane) [N]
    F_cor     = mass_vane*aGcor_j;                % Corolis force (perpendicular to vane) [N]
    C_izG     = moment_inertia_vane * DomegaG;    % inertial torque [N m]
    % =====================================================================
    % =====================================================================
    Fw_i      =  rot_dir*mass_vane*g_const*sin(theta_kin)*cos(zeta)-mass_vane*g_const*cos(theta_kin)*sin(zeta);  % i-comp. weigth force [N]
    Fw_j      = -rot_dir*mass_vane*g_const*cos(theta_kin)*cos(zeta)-mass_vane*g_const*sin(theta_kin)*sin(zeta);  % j-comp. weigth force [N]
    % =====================================================================
    % =====================================================================


    %% PRELIMINARY SOLUTION - Cramer analytical method %%
    %(taken from: Analysis of Mechanical Friction in rotany Vane machines(1972)- Edwards;McDonald)
    
    % Definition of coefficients and matrix A elements for configuration 1
    f_cv = f_c*sign(vAGi);                        %  array of slot friction coefficient
    a11 = -cos(UPSILON)- f_t*sin(UPSILON);
    a12 = -f_cv;
    a13 = -f_cv;
    a21 = -sin(UPSILON)+ f_t*cos(UPSILON);
    a22 = -1;
    a23 = +1;
    a31 = +sin(UPSILON).*AP_i - f_t*cos(UPSILON).*AP_i - cos(UPSILON).*AP_j - f_t*sin(UPSILON).*AP_j;
    a32 = AM_i(1,:) - 0.5*s*f_cv;
    a33 = -AV_i + 0.5*s*f_cv;

    % Definition of vector of known elements for configuration 1
    b1 =  (F_iz_i - F_centr_i + Fw_i) - Fp_cava + Fp_1_tip.*cos(UPSILON_tip_1) + Fp_2_tip.*cos(UPSILON_tip_2);
    b2 =  -(sgn*(F_cor  + F_centr_j) - Fw_j) - Fp_2 - Fp_2_cava + Fp_1 + Fp_1_cava + Fp_1_tip.*sin(UPSILON_tip_1) + Fp_2_tip.*sin(UPSILON_tip_2);   % (SEE NOTE #4)
    b3 =  (sgn*(F_cor  + F_centr_j) - Fw_j).*AG + (F_centr_i - F_iz_i - Fw_i)*b_gx - Fp_1.*AX_i - Fp_1_cava.*Ax_i + Fp_2.*AY_i + Fp_2_cava.*Ay_i +...
          Fp_1_tip.*(-sin(UPSILON_tip_1).*AR_i + cos(UPSILON_tip_1).*AR_j) + Fp_2_tip.*(-sin(UPSILON_tip_2).*AS_i + cos(UPSILON_tip_2).*AS_j) + C_izG;
    
    % plot di verifica
    % Vettore posizioni angolari
    %theta_pt =(1:Npt);
    %theta_deg = theta_pt.*360/Npt;

    %figure, hold on,grid on, box on;
    %plot(theta_deg,UPSILON,'g',theta_deg,UPSILON_tip_1,'b',theta_deg,UPSILON_tip_2,'r');
    %xlim([0 360]);
    %title('Verifica');
    %xlabel('theta [deg]');
    %ylabel('[-]');
    %legend('UPSILON','1 (low)','2 (high)');
    %fine plot di verifica

    % Calculation of determinants of matrixes A, A_Fm and A_Fv
    detA    = a11.*a22.*a33 + a12.*a23.*a31 + a13.*a21.*a32 - a31.*a22.*a13 - a32.*a23.*a11 - a33.*a21.*a12;
    detA_Fm = a11.*b2.*a33 + b1.*a23.*a31 + a13.*a21.*b3 - a31.*b2.*a13 - b3.*a23.*a11 - a33.*a21.*b1;
    detA_Fv = a11.*a22.*b3 + a12.*b2.*a31 + b1.*a21.*a32 - a31.*a22.*b1 - a32.*b2.*a11 - b3.*a21.*a12; 

    % Calculation of unknowns F_mAp e F_vAp for configuration 1 (preliminary solutions)
    F_mAp = detA_Fm./detA;
    F_vAp = detA_Fv./detA;
    F_vAp(F_vAp==0) = 1e-12;  % top avoid null value, otherwise sign() is equal to zero!
    F_vAp(F_mAp==0) = 1e-12;

    %% CONFIGURATIONS IDENTIFICATION - analytical method %%
    s1   = NaN(1,Npt);      % sign of F_m
    s2   = NaN(1,Npt);      % sign of F_v
    conf = NaN(1,Npt);      % configuration array

    % I configuration
    cond1       = sign(F_mAp) == 1 & sign(F_vAp) == 1;
    s1(cond1)   = 1;
    s2(cond1)   = 1;
    conf(cond1) = 1;
    % II configuration
    cond2       = sign(F_mAp) == -1 & sign(F_vAp) == -1; 
    s1(cond2)   = -1;
    s2(cond2)   = -1;
    conf(cond2) =  2;
    % III configuration
    cond3       = sign(F_mAp) == -1 & sign(F_vAp) == 1; 
    s1(cond3)   = -1;
    s2(cond3)   =  1;
    conf(cond3) =  3;
    % IV configuration
    cond4       = sign(F_mAp) == 1 & sign(F_vAp) == -1;
    s1(cond4)   =  1;
    s2(cond4)   = -1;
    conf(cond4) =  4;
    clear F_mAp F_vAp

    % angle at which there is a change in configuration
    idx     = diff(conf)~=0;
    th_conf = theta_vane(idx);
    
    %% REAL SOLUTION - Cramer analytical method %%
    % correction of  AM_i vector depending on configuration
    AM_2 = AM_i(2,:);
    AM_3 = AM_i(3,:);
    AM_4 = AM_i(4,:);

    AM_i_new = AM_i(1,:);
    AM_i_new(conf == 2) = AM_2(conf == 2);
    AM_i_new(conf == 3) = AM_3(conf == 3);
    AM_i_new(conf == 4) = AM_4(conf == 4);
    clear AM_2 AM_3 AM_4 

    % Correction of matrix A elements depending on configuration
    a12 = - s1.*f_cv;
    a13 = - s2.*f_cv;
    a32 = AM_i_new - 0.5*s*s1.*f_cv.*s1;  %(SEE NOTE #5)
    a33 = -AV_i + 0.5*s*s2.*f_cv.*s2;     %(SEE NOTE #5)

    % Calculation of determinants of corrected matrixes A, A_Ft, A_Fm and A_Fv
    detA    = a11.*a22.*a33 + a12.*a23.*a31 + a13.*a21.*a32 - a31.*a22.*a13 - a32.*a23.*a11 - a33.*a21.*a12;
    detA_Ft = b1.*a22.*a33  + a12.*a23.*b3  + a13.*b2.*a32  - b3.*a22.*a13  - a32.*a23.*b1  - a33.*b2.*a12;
    detA_Fm = a11.*b2.*a33  + b1.*a23.*a31  + a13.*a21.*b3  - a31.*b2.*a13  - b3.*a23.*a11  - a33.*a21.*b1;
    detA_Fv = a11.*a22.*b3  + a12.*b2.*a31  + b1.*a21.*a32  - a31.*a22.*b1  - a32.*b2.*a11  - b3.*a21.*a12; 

    % Calculation of corrected unknowns F_tA, F_mAp and F_vAp [N]
    F_tA = detA_Ft./detA; % Reazione vincolare al tip
    F_mA = detA_Fm./detA; % Reazione vincolare a monte
    F_vA = detA_Fv./detA; % Reazione vincolare a valle 

    % Friction forces are computed [N]
    T_tA = f_t*F_tA;       % attrito generato da F_tA
    T_mA = s1.*f_cv.*F_mA; % attrito generato da F_mA
    T_vA = s2.*f_cv.*F_vA; % attrito generato da F_vA

    % torque due to axial forces
    Mf_ax = s2.*T_vA*0.5*s - s1.*T_mA*0.5*s - sgn*F_tA.*cos(UPSILON).*AP_j - sgn*T_tA.*sin(UPSILON).*AP_j + (F_centr_i - F_iz_i - Fw_i)*b_gx + Fp_1_tip.*sgn.*cos(UPSILON_tip_1).*AR_j + Fp_2_tip.*sgn.*cos(UPSILON_tip_2).*AS_j;
    
    clear cond1 cond2 cond3 cond4
    clear a11 a12 a13 a21 a22 a23 a31 a32 a33 detA detA_Ft det_Fm det_Fv

    %% PRELIMINARY SOLUTION - numerical method %%
    if fDBG == 1    
        % Preallocation
        A     = NaN(3,3,Npt); % coefficent matrix
        b     = NaN(3,1,Npt); % known elemetns vector
        F     = NaN(3,1,Npt); % unknown forces vector

        % Preliminary system definition (unknowns: F_t, F_m, F_v)
        A(1,1,:) = -cos(UPSILON)- f_t*sin(UPSILON);
        A(1,2,:) = -f_cv;
        A(1,3,:) = -f_cv;
        A(2,1,:) = -sin(UPSILON)+ f_t*cos(UPSILON);
        A(2,2,:) = -1;
        A(2,3,:) = +1;
        A(3,1,:) = +sin(UPSILON).*AP_i - f_t*cos(UPSILON).*AP_i - cos(UPSILON).*AP_j - f_t*sin(UPSILON).*AP_j;
        A(3,2,:) = AM_i(1,:) - 0.5*s*f_cv;
        A(3,3,:) = -AV_i + 0.5*s*f_cv;
        b(1,1,:) = (F_iz_i - F_centr_i + Fw_i) - Fp_cava + Fp_1_tip.*cos(UPSILON_tip_1) + Fp_2_tip.*cos(UPSILON_tip_2);
        b(2,1,:) = -(sgn*(F_cor  + F_centr_j) - Fw_j) - Fp_2 - Fp_2_cava + Fp_1 + Fp_1_cava + Fp_1_tip.*sin(UPSILON_tip_1) + Fp_2_tip.*sin(UPSILON_tip_2); % (SEE NOTE #4)
        b(3,1,:) = (sgn*(F_cor  + F_centr_j) - Fw_j).*AG + (F_centr_i - F_iz_i - Fw_i)*b_gx - Fp_1.*AX_i - Fp_1_cava.*Ax_i + Fp_2.*AY_i + Fp_2_cava.*Ay_i +...
                   Fp_1_tip.*(-sin(UPSILON_tip_1).*AR_i + cos(UPSILON_tip_1).*AR_j) + Fp_2_tip.*(-sin(UPSILON_tip_2).*AS_i + cos(UPSILON_tip_2).*AS_j) + C_izG;

        % Preliminary system solution & solution extraction
        for i=1:Npt
            F(:,:,i)=A(:,:,i)\b(:,:,i);
        end, clear i
        F_mNp(1:Npt) = F(2,1,:);
        F_vNp(1:Npt) = F(3,1,:);

        %% CONFIGURATIONS IDENTIFICATION - numerical method %%
        % identificatioon of configuration and correction of matrix A
        % II configurazione
        cond2 = sign(F_mNp) == -1 & sign(F_vNp) == -1;
        AM_i_conf2   = AM_i(2,:);
        A(1,2,cond2) = f_cv(cond2);
        A(1,3,cond2) = f_cv(cond2);
        A(3,2,cond2) = AM_i_conf2(cond2) + 0.5*s*f_cv(cond2)*(-1);
        A(3,3,cond2) = -AV_i(cond2) - 0.5*s*f_cv(cond2)*(-1);
        % III configurazione
        cond3 = sign(F_mNp) == -1 & sign(F_vNp) == 1; 
        AM_i_conf3   = AM_i(3,:);
        A(1,2,cond3) = f_cv(cond3);
        A(1,3,cond3) = -f_cv(cond3);
        A(3,2,cond3) = AM_i_conf3(cond3) + 0.5*s*f_cv(cond3)*(-1);
        A(3,3,cond3) = -AV_i(cond3) + 0.5*s*f_cv(cond3);
        % IV configurazione
        cond4 = sign(F_mNp) == 1 & sign(F_vNp) == -1;
        AM_i_conf4   = AM_i(4,:);
        A(1,2,cond4) = -f_cv(cond4);
        A(1,3,cond4) = f_cv(cond4);
        A(3,2,cond4) = AM_i_conf4(cond4) - 0.5*s*f_cv(cond4);
        A(3,3,cond4) = -AV_i(cond4) - 0.5*s*f_cv(cond4)*(-1);

        %% REAL FORCES CALCULATIONS - numerical method %%
        % system solve
        for i=1:Npt
             F(:,:,i)=A(:,:,i)\b(:,:,i);
        end, clear i

        % Real solution calculation
        F_tN(1:Npt) = F(1,1,:);  % tip reaction force [N]
        F_mN(1:Npt) = F(2,1,:);  % reaction force at rotor-vane contact point [N]
        F_vN(1:Npt) = F(3,1,:);  % reaction force at vane bottom [N]

    end, clear A b F F_mNp F_vNp cond1 cond2 cond3 cond4 method

    %% CHECKS %%
    % Check on dynamic equilibrium of axial and perpendicular forces, and of moment respect to point A
    Eql_N  = -T_mA -T_vA - F_tA.*cos(UPSILON) - T_tA.*sin(UPSILON) - b1;
    Eql_T  = F_vA-F_mA - F_tA.*sin(UPSILON) +T_tA.*cos(UPSILON) - b2;
    Eql_Mf = -AV_i.*F_vA + s2.*T_vA*0.5*s + F_mA.*AM_i_new - s1.*T_mA*0.5*s + F_tA.*sin(UPSILON).*AP_i - sgn*F_tA.*cos(UPSILON).*AP_j - T_tA.*cos(UPSILON).*AP_i - sgn*T_tA.*sin(UPSILON).*AP_j - b3; 
    
    chk_N  = abs(Eql_N)> toll_t;
    chk_T  = abs(Eql_T)> toll_t; 
    chk_Mf = abs(Eql_Mf) > toll_t;
    
    if sum(chk_N)>=1
        maxN = max(abs([T_mA, T_vA,F_tA.*cos(UPSILON),T_tA.*sin(UPSILON), b1]));
        warning('S4_VaneDynamics_C:Equilibrium','Vane axial forces not in equilibrium. Estimated relative error: %2.6f %%',max(chk_N)/maxN*100);
        SX_Logfile ('w',{lastwarn});
    end
    if sum(chk_T)>=1
        maxT = max(abs([F_vA, F_mA F_tA.*sin(UPSILON), T_tA.*cos(UPSILON), b2]));
        warning('S4_VaneDynamics_C:Equilibrium','Vane perpendicular forces not in equilibrium. Estimated relative error: %2.6f %%',max(chk_T)/maxT*100);
        SX_Logfile ('w',{lastwarn});
    end
    if sum(chk_Mf)>=1
        maxMf = max(abs([AV_i.*F_vA, s2.*T_vA*0.5*s, F_mA.*AM_i_new, s1.*T_mA*0.5*s, F_tA.*sin(UPSILON).*AP_i, sgn*F_tA.*cos(UPSILON).*AP_j, T_tA.*cos(UPSILON).*AP_i, sgn*T_tA.*sin(UPSILON).*AP_j, b3 ]));
        warning('S4_VaneDynamics_C:Equilibrium','Vane bending moment not in equilibrium. Estimated relative error: %2.6f %%',max(chk_Mf)/maxMf*100);
        SX_Logfile ('w',{lastwarn});
    end, clear chk_N chk_T chk_Mf Eql_N Eql_T Eql_Mf
    
    
    % Check on vane-stator contact
    chck_tip = F_tA<0;
    if sum(chck_tip)>=1
        warning('S4_VaneDynamics_C:Contact','No contact between vane tip and stator')
        SX_Logfile ('w',{lastwarn});
    end
    chck_tip = ~chck_tip;
    
    % Comparison between analitical and nunumerical model
    if fDBG == 1    
        errF_t = abs((F_tN-F_tA)./F_tA);
        errF_m = abs((F_mN-F_mA)./F_mA);
        errF_v = abs((F_vN-F_vA)./F_vA);
        
        chck1 = errF_t>toll_t;
        chck2 = errF_m>toll_t;
        chck3 = errF_v>toll_t;
        chck4 = isnan(conf);
        
        if sum(chck1)>=1
            warning('S4_VaneDynamics_C:forces','Analytical and numerical model for tip forces do not match. Max relative error: %2.4f %%', max(errF_t)*100);
            SX_Logfile ('d',{lastwarn});
        end
        if sum(chck2)>=1
            warning('S4_VaneDynamics_C:forces','Analytical and numerical model for middle-vane forces do not match. Max relative error: %2.4f %%', max(errF_m)*100);
            SX_Logfile ('d',{lastwarn});
        end
        if sum(chck3)>=1
            warning('S4_VaneDynamics_C:forces','Analytical and numerical model for vane-bottom forces do not match. Max relative error: %2.4f %%', max(errF_v)*100);
            SX_Logfile ('d',{lastwarn});
        end
        if sum(chck4)>=1
            warning('S4_VaneDynamics_C:forces','Invalid vane configuration')
            SX_Logfile ('d',{lastwarn});
        end, clear errF_t errF_m errF_v chck1 chck2 chck3 chck4 F_tN F_mN F_vN
    end       

end
