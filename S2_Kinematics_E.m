function [mass_vane,Rth,xG,w,beta,psi,vGr_A,aGr_iz_A,aGr_centr,aGt_cor,omega,sgn] = S2_Kinematics_E(d,e,s,L,rho_vane,theta_vane,rpm,l,toll_d,Npt,fDBG)
% Calcola le grandezze geometriche, la velocità e le accelerazioni del
% baricentro per statore ellittico
% INPUT
% d [m]            : rotor diameter
% e [-]            : ellipse eccentricity
% s [m]            : vane thickness
% L [m]            : rotor lenght
% rho_vane [kg/m3] : vane density 
% theta_vane[rad]  : vane angular positions between (0:2pi)
% rpm [rpm]        : angular velocity
% l [m]            : vane heigth
% toll_d [-]       : numeric tolerance
% Npt [-]          : discretization points
% fDBG [bool]      : denug flag
%
% OUTPUT
% mass_vane [kg]    : mass of one vane
% Rth [m]           : stator radius seen by rotor center O
% xG [m]            : array of center of mass position
% w [m]             : vane excursion
% beta [rad]        : pressure angle at vane tip
% psi [m]           : BOC angle (it is zero for radial vanes)
% vGr_A [m/s]       : array of radial component of center of mass velocity (analytical)
% aGr_iz_A [m/s^2]  : array of radial component of center of mass inertial acceleration (analytical)
% aGr_centr [m/s^2] : array of radial component of center of mass centripetal acceleration (analytical)
% aGt_cor [m/s^2]   : array of tangential component of center of mass Coriolis acceleration (analytical)
% omega [rad/s]     : angular velocity
% sgn [-]           : vane tilt direction
%
% HP: - zero thickness vanes
%     - radial vanes

    %% GEOMETRY CALCULATIONS %%
    b         = d/2;                                           % short axes [m]
    a         = b/sqrt(1-e^2);                                 % long axes [m]
    Rth       = a*sqrt(1-e^2)./sqrt(1-(e*sin(theta_vane)).^2); % stator radius seen by rotor center O [m]
    omega     = 2*pi*rpm/60;                                   % angular velocity [rad/s]
    xG        = Rth-l/2;                                       % array of center of mass position [m]
    w         = Rth-d/2;                                       % vane excursion [m]
    psi       = 0;                                             % BOC angle (it is zero for radial vanes)
    sgn       = 1;                                             % tilted vane discriminant factor (it is one for radial vanes) 
    mass_vane = rho_vane*s*L*l;                                % vane mass è [kg]
    
    
    % Pressure angle at vane tip
    % NB: eliptical vane are 1D -> no tip radius is considered
    % (taken from scientific article);it just find the perpendicular direction to stator
    beta = atan(0.5*e^2*sin(2*theta_vane)./(1-(e*sin(theta_vane)).^2));

    %% KINEMATICS COMPUTATIONS - analytical method %%
    % definition of temporary factors
    dt   = 60/rpm/Npt/2;   % discretized time step [s]
    Ft1  = a*sqrt(1-e^2)*e^2;
    Ft2  = (1-(e*sin(theta_vane)).^2).^(3/2);
    t2_p = 2*theta_vane; 

    % analytical velocity [m/s]
    vGr_A      = omega*(Ft1/2*sin(2*theta_vane)./Ft2); % radial velocity
    vGr_A(1)   = 0;                                    % in tangency, radial velocity is null
    vGr_A(end) = 0;
    % vGt_A    = xG*omega;                             % tangential velocity

    % analytical acceleration [m/s2]
    aGr_iz_A     = omega^2*((Ft1./Ft2).*((3/4*e^2*(sin(t2_p)).^2)./(Ft2.^(2/3))+cos(t2_p)));  % radial acceleration
    aGr_centr = xG*(omega^2);                                                              % centrifugal acceleration (radial)
    aGt_cor   = 2*vGr_A*omega;                                                             % Coriolis acceleration (tangential)
    clear Ft1 Ft2 t2_p

    %% KINEMATIC COMPUTATIONS - numerical method %%
    if fDBG == 1
        i   = 2:Npt-1;
        % creation of numerical derivative using central difference method (II order method)
        vGr_N = ((xG(i+1)-xG(i-1))/(2*dt));     % numerical velocity of G [m/s]
        aGr_N = (vGr_A(i+1)-vGr_A(i-1))/(2*dt); % numerical acceleration of G [m/s2] 
    end

        %% CHECKS %%
        if fDBG == 1
            err_v = abs((vGr_N-vGr_A(2:end-1))./vGr_A(2:end-1));       % velocity relative error
            err_a = abs((aGr_N-aGr_iz_A(2:end-1))./aGr_iz_A(2:end-1)); % acceleration relative error
            
            % infinite and unitary values are canceled
            err_v(isinf(err_v)) = [];
            err_v(err_v == 1)   = [];
            err_a(isinf(err_a)) = [];
            err_a(err_a == 1)   = [];
            
            % warning display
            chck1=err_v > toll_d;
            chck2=err_a > toll_d;

            if sum(chck1) >= 1
                warning('S2_Kinematics_E:Velocity','Velocity: numerical and analytical model do not match. Max relative error: %2.4f %% ', max(err_v)*100)
                SX_Logfile ('d',{lastwarn});
            end
            if sum(chck2) >= 1
                warning('S2_Kinematics_E:accelerattion','Acceleration: numerical and analytical model do not match. Max relative error: %2.4f %% ', max(err_a)*100)
                SX_Logfile ('d',{lastwarn});
            end        
            clear err_v chck1 err_a chck2 i
        end

end