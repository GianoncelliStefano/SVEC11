function [AG,b_gx,xG,rv,psi,alpha,vAGi_A,aGi_A,aGcentr_i,aGcentr_j,aGcor_j,aIVcentr_i,omega,DomegaG_A]...
    = S2_Kinematics_C(D,d,b,Y_F,Xg_vane,Yg_vane,R_v,theta_vane,rpm,r_tip,s,l,th_tilt,sgn,AV_i,AI,Npt,toll_t,toll_d,fDBG)
% Calcola le grandezze geometriche e risolve la cinematica del baricentro
%
% INPUT
% D [m]             : stator diameter
% d [m]             : rotor diameter
% b [%]             : tip radius offset
% Y_F [m]           : y-coord of vane tip center of curvature
% Xg_vane [m]       : x-coord of vane center of mass
% Yg_vane [m]       : y-coord of vane center of mass
% R_v [m]           : perpendicular to vane passing from stator center
% theta_vane [rad]  : array of discretized angular position of trailing vane [-gamma:2pi]
% rpm [giri/min]    : rotational speed
% r_tip [m]         : tip radius
% s [m]             : vane thickness
% th_tilt [rad]     : vane tilt angle
% sgn [-]           : vane tilt angle factor
% AV_i, AI [m]      : geometrical quantities (see D Zanchi Thesis)
% Npt [-]           : number of discretization points
% toll_t [-]        : theoric tolerance
% toll_d [-]        : numeric tolerance
% fDBG              : developer control flag
%
% OUTPUT
% AG [m]             : distance of vane center of mass respect to point A
% b_gx [m]           : offset of vane center of mass along vane axis
% xG [m]             : position of vane center of mass
% rv [m]             : OA segment
% psi [rad]          : BOC angle
% alpha [rad]        : angular position of point A
% vAGi_A [m/s]       : axial component of vane center of mass velocity vector (analitical)
% aGi_A [m/s2]       : axial component of vane center of mass inertial acceleration (analitical)
% aGcentr_i [m/s2]   : axial component of vane center of mass centripetal acceleration (axial to vane)
% aGcentr_j [m/s2]   : perpendicular component of vane center of mass centripetal acceleration (perpendicular to vane)
% aGcor_j [m/s2]     : perpendicular component of vane center of mass Coriolis acceleration (perpendicular to vane)
% aIVcentr_i [m/s2]  : matrix that has on the first row the vector of the axial centripetal acceleration of point I (vane tip)
%                       and on the second row the the vector of ther axial centripetal acceleration of point V (vane bottom)
% omega [rad/s]      : angular speed
% DomegaG_A [rad/s2] : angular accelleration of vane center of mass (analitical)

% NOTES: 
% the kinematics quantities are referred to the TRAILING vane
 
    %% DECLARATIONS %%
    r     = d/2;           % rotor radius [m]
    R     = D/2;           % stator radius [m]
    e     = R-r;           % eccentricity of machine [m]
    omega = 2*pi*rpm/60;   % angular velocity [rad/s]

    %% GEOMETRIC DEFINITIONS %%
    rv   = abs(r*sin(th_tilt));        % OA segment  [m]
    eta  = sgn*(pi/2- abs(th_tilt));   % AOB angle [rad]
    b_gx  = 0.5*s - Xg_vane;           % vane center of mass position respect to vane mid plane [m]
    b_gy = Y_F - Yg_vane;              % vane center of mass position respect to tip radius center [m]

    alpha    = theta_vane - eta;                                % [rad]
    AH       = sqrt((R-r_tip)^2 - R_v.^2) + sgn*e*sin(alpha);   % [m]
    AC       = sgn*e*sin(alpha)+sqrt(R^2-(rv+e*cos(alpha)).^2); % [m] old geometry model 2D vane
    AG       = AH - b_gy;                                       % [m]
    OC       = sqrt(AC.^2+ (rv)^2);                             % [m]
    AIV      = NaN(2,length(theta_vane));                       %      preallocation
    AIV(1,:) = AI;                                              % [m]  positions of vane tip
    AIV(2,:) = AV_i;                                            % [m]  position of vane bottom
    xG       = sqrt(AG.^2+(rv-abs(b_gx))^2);                    % [m]  vane center of mass distance respect to rotor center - correspond to OG
    nu       = asin((rv-abs(b_gx))./xG);                        % [rad] OGA angle
    phi      = asin(AC./OC);                                    % [rad] AOC angle
    psi      = phi-abs(eta);                                    % [rad] BOC angle
    betaN = asin((e*sin(pi-theta_vane-sgn*psi))/R);             % [rad] OCS abgle
    delta    = asin(AG./xG)- sgn* eta;                          % [rad] BOG angle

    %% VANE KINEMATICS: ANALYTIC CALCULATIONS %%
    % Definition of temporary factor useful for derivative calculation
    theta_kin   = theta_vane + th_tilt;
    I1          = e*cos(theta_kin)*omega;
    I2          = (e+r*cos(theta_vane)).*sin(theta_kin) - r*sin(theta_vane).*cos(theta_kin) - b;
    f           = I1.*I2;
    g           = sqrt((R-r_tip)^2 - I2.^2);
    f_primo     = I1.^2 + I2.*(-e*sin(theta_kin).*omega^2);
    g_primo     = -(I1.*I2)./g;
    der_fg = ((f_primo.*g-f.*g_primo)./(g.^2));

    % Calculation of center of mass velocity - analytical method [m/s] o [rad/s]
    vAGi_A   = (sgn*e*cos(alpha)*omega-(f./g));         % relative velocity of G respect to A - i component
    vAGj_A   = sgn*AG*omega;                            % relative velocity respect of G to A - j component
    vAi_A    = sgn*(rv-abs(b_gx))*omega;                % velocity of A respect to rotor center - i component
    vGi_A    = vAGi_A + vAi_A;                          % velocity of G respect to rotor center - i component
    vG_t     = sgn.*(vGi_A.*sin(nu) + vAGj_A.*cos(nu)); % tangential velocity of G (perpendicular to vector xG)
    vG_r     = vGi_A.*cos(nu) - vAGj_A.*sin(nu);        % radial velocity of G (parallel to vector xG)
    omegaG_A = vG_t./xG;                                % angular velocity of G

    % Calculation of center of mass acceleration - analytical method [m/s2]
    aGi_A     = (-sgn*e*sin(alpha).*omega^2 - der_fg); % acc. of G - i component
    aGcentr_i = AG.*omega.^2;                          % centripetal acc. of G -  i component
    aGcentr_j = (rv-abs(b_gx))*omega.^2;               % centripetal acc. of G -  j component
    aGcor_j   = sgn.*2*vAGi_A.*omega;                  % coriolis acc. of G -  j component
    aGcor_t   = 2.*omegaG_A.*vG_r;                     % coriolis acc. of G -  tangential component
    aGri_t    = sgn.*aGi_A.*sin(nu);                   % acc. of G -  tangential component
    aGj_cor   = sgn.*aGcor_j.*cos(nu);                 %
    aG_t      = aGj_cor + aGri_t - aGcor_t;            % acc. tangenziale di G

    % Calculation of G angular acceleration [rad/s2]
    if th_tilt == 0
        DomegaG_A = zeros(1,length(theta_vane)); 
    else
        DomegaG_A = aG_t./xG;
    end
    % Acceleration of vane bottom and tip points
    aIVcentr_i  = AIV.*omega.^2; % [m/s^2]
    clear AB eta b_gy AH OC ACD theta_kin I1 I2 f g f_primo g_primo der_ratiofg vAGj_A vAi_A vG_t vG_r aGcor_t aGri_t aGj_cor aG_t

    %% CHECKS %%
    check_1=sum(AG<0);
    if check_1>0
        warning('S2_Geometry_C:accelerazione','AG is negative for some angular positions')
        SX_Logfile ('w',{lastwarn});
    end
    clear check_1

    %% VANE KINEMATICS: NUMERIC SOLUTIONS %%
    if fDBG == 1
        dt = 60/rpm/Npt;          % time step (considering 360°) [s]
        i= 2:Npt-1;
        % creation of numerical derivative using central difference method (II order method)
        vAGi_N    = ((AG(1,i+1) - AG(1,i-1))/(2*dt));           % numerical velocity of G [m/s]
        aGi_N     = (vGi_A(i+1) - vGi_A(i-1))/(2*dt);           % numerical acceleration of G [m/s2]
        omegaG_N  = omega+sgn*((delta(i+1)-delta(i-1))/(2*dt)); % numerical angular velocity of G [rad/s]
        DomegaG_N = (omegaG_A(i+1)-omegaG_A(i-1))/(2*dt);       % numerical angular acc. of G [rad/s2]

        % creation of error vectors
        err_v  = abs((vAGi_N-vAGi_A(2:end-1))./vAGi_A(2:end-1));           % relative error on G velocity
        err_a  = abs((aGi_N-aGi_A(2:end-1))./aGi_A(2:end-1));              % relative error on G acceleration
        err_dv = abs((omegaG_N-omegaG_A(2:end-1))./omegaG_A(2:end-1));     % relative error on G angular velocity 
        err_da = abs((DomegaG_N-DomegaG_A(2:end-1))./DomegaG_A(2:end-1));  % relative error on G angular acceleration
        
        % infinite and values are canceled
        err_v(isinf(err_v))   = [];
        err_v(err_v == 1)     = [];
        err_a(isinf(err_a))   = [];
        err_a(err_a == 1)     = [];
        err_dv(isinf(err_dv)) = [];
        err_dv(err_dv == 1)   = [];
        err_da(isinf(err_da)) = [];
        err_da(err_da == 1)   = [];
        
        % warning display
        chck1=err_v > toll_d;
        chck2=err_a > toll_d;
        chck3=err_dv > toll_d;
        chck4=err_da > toll_d;
        
        if sum(chck1) >= 1
            warning('S2_Kinematics_C:Velocity','Velocity: numerical and analytical model do not match. Max relative error: %2.4f %% ', max(err_v)*100)
            SX_Logfile ('d',{lastwarn});
        end
        if sum(chck2) >= 1
            warning('S2_Kinematics_C:Acceleration','Acceleration: numerical and analytical model do not match. Max relative error: %2.4f %% ', max(err_a)*100)
            SX_Logfile ('d',{lastwarn});
        end
        if sum(chck3) >= 1
            warning('S2_Kinematics_C:Vel_angular','Angular velocity: numerical and analytical model do not match. Max relative error: %2.4f %% ', max(err_dv)*100)
            SX_Logfile ('d',{lastwarn});
        end
        if sum(chck4) >= 1
            warning('S2_Kinematics_C:Acc_angular','Angular acc.: numerical and analytical model do not match. Max relative error: %2.4f %% ', max(err_da)*100)
            SX_Logfile ('d',{lastwarn});
        end
        clear i chck1 chck2 chck3 chck4 err_v err_a err_da err_dv delta aPi_N vAPi_N omegaG_A omegaG_N DomegaG_N
    end

    %% COMPARISON BETWEEN NEW AND OLD MODELS WHEN b,r_tip and s ARE BROUGHT TO ZERO %%
    if fDBG == 1
        if b == 0 && r_tip<toll_t && s<toll_t
            AG_Old = AC - l/2;             % [m]
            xG_Old = sqrt(AG_Old.^2+rv^2); % [m]
            nu_Old = asin(rv./xG_Old);     % [rad]

            % Definition of temporary factor useful for derivative calculation
            Ft1 = sqrt(R^2-(rv+e*cos(alpha)).^2);  % [m]
            Ft2 = e*sin(alpha).*(e*cos(alpha)+rv); % [m^2]
            Ft3 = Ft2./Ft1 +sgn*e*cos(alpha);      % [m]
            Ft4 = e*cos(alpha).*(e*cos(alpha)+rv); % [m^2]
            Ft5 = -sgn*e*sin(alpha);               % [m]
            Ft6 = (Ft2.^2)./(Ft1.^3);              % [m]
            
            % Calculation of center of mass velocity - analytical method [m/s] o [rad/s]
            vAGi_A_Old = Ft3*omega;                                                % relative velocity of G respect to A - i component
            vAGj_A_Old = sgn*AG_Old*omega;                                         % relative velocity of G respect to A - j component
            vAi_A_Old  = sgn*rv*omega;                                             % velocity of A respect to rotor center - i component
            vGi_A_Old  = vAGi_A_Old + vAi_A_Old;                                   % velocity of G respect to rotor center - i component
            vG_t_Old   = sgn.*(vGi_A_Old .*sin(nu_Old) + vAGj_A_Old.*cos(nu_Old)); % tangential velocity of G (perpendicular to vector xG)
            vG_r_Old   = vGi_A_Old .*cos(nu_Old) - vAGj_A_Old.*sin(nu_Old);        % radial velocity of G (parallel to vector xG)
            omegaG_A_Old = vG_t_Old./xG_Old;                                       % angular velocity of G
            
            % Calculation of center of mass acceleration - analytical method [m/s2]
            aGi_A_Old     = ( Ft5 + (Ft4 - Ft5.^2)./Ft1 - Ft6 )*omega^2;  % acc. of G - i component
            aGcentr_i_Old = AG_Old.*omega.^2;                             % centripetal acc. of G -  i component
            aGcentr_j_Old = rv*omega.^2;                                  % centripetal acc. of G -  j component
            aGcor_j_Old   = sgn.*2*vAGi_A_Old.*omega;                     % coriolis acc. of G -  j component
            aGcor_t_Old   = 2.*omegaG_A_Old.*vG_r_Old;                    % coriolis acc. of G -  tangential component
            aGri_t_Old    = sgn.*aGi_A_Old.*sin(nu_Old);                  % acc. of G -  tangential component
            aGj_cor_Old   = sgn.*aGcor_j_Old.*cos(nu_Old);                % 
           % Calculation of G angular acceleration [rad/s2]
            if th_tilt == 0
                aG_t_Old = 0;
            else
                aG_t_Old = aGj_cor_Old + aGri_t_Old - aGcor_t_Old;
            end

            % Angular acceleration [rad/s2]
            DomegaG_A_Old = aG_t_Old./xG_Old;

            % creation of error vectors
            err_1 = abs((aGi_A-aGi_A_Old)./min(aGi_A,aGi_A_Old));                 % relative error on G acc.
            err_2 = abs((aGcentr_i-aGcentr_i_Old)./min(aGcentr_i,aGcentr_i_Old)); % relative error on G centrifugal acc.
            err_3 = abs((aGcentr_j-aGcentr_j_Old)./min(aGcentr_j,aGcentr_j_Old)); % relative error on G centrifugal acc.
            err_4 = abs((DomegaG_A-DomegaG_A_Old)./min(DomegaG_A,DomegaG_A_Old)); % relative error on G angular acc.
            err_5 = abs((vGi_A-vGi_A_Old)./min(vGi_A,vGi_A_Old));                 % relative error on G velocity
            chck1 = err_1 > toll_t;
            chck2 = err_2 > toll_t;
            chck3 = err_3 > toll_t;
            chck4 = err_4 > toll_t;
            chck5 = err_5 > toll_t;
            % Unitary values are canceled
            chck1(err_1==1) = [];
            chck2(err_2==1) = [];
            chck3(err_3==1) = [];
            chck4(err_4==1) = [];
            chck5(err_5==1) = [];

            % Display Warnings
            if sum(chck1) >=1
                warning('S2_Kinematics_C:acceleration','Accelerations of 1D and 2D models does not match')
                SX_Logfile ('d',{lastwarn});
            end

            if sum(chck2) >=1
                warning('S2_Kinematics_C:acceleration','Centrifugal acceleration of 1D and 2D models does not match')
                SX_Logfile ('d',{lastwarn});
            end

            if sum(chck3) >=1
                warning('S2_Kinematics_C:acceleration','Centrifugal acceleration of 1D and 2D models does not match')
                SX_Logfile ('d',{lastwarn});
            end

            if sum(chck4) >=1
                warning('S2_Kinematics_C:acceleration','Angular acceleration of 1D and 2D models does not match')
                SX_Logfile ('d',{lastwarn});
            end

            if sum(chck5) >=1
                warning('S2_Kinematics_C:velocity','Velocity of 1D and 2D models does not match')
                SX_Logfile ('d',{lastwarn});
            end

        end
    end
end
 

