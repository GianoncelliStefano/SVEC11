function V_cell = S2_Volume_C(D,d,L,theta,Gamma,sgn,Npt_cell,pos_tang,th_tilt,toll_d,fDBG,Vol_Excluded,s)
% this function calculates the volumes with theta that goes from 0-gamma to 2pi
% INPUT
% D             [m]     : stator diameter
% d             [m]     : rotor diameter
% L             [m]     : compressor axial length
% theta         [rad]   : vane angular positions between (-gamma:2pi)
% Gamma         [rad]   : angular distance between two vanes
% sgn           [-]     : sign of vane tilt angle
% Npt_cell      [-]     :  number of discretization points of a cell (space between 2 vanes)
% pos_tang      [-]     : tangency position
% th_tilt       [rad]   : vanes inclination angle
% toll_d        [-]     : numeric tolerance
% fDBG                  : debug flag
% Vol_Excluded          : (3,:) matrix of excluded volume due to 3D vane
%                       : 1° row is first half vane volume
%                       : 2° row is second half vane volume
%                       : 3° row is volume between vane tip and stator
%
% OUTPUT
% V_cell        [m^3]   : cell volumes array of a complete vane revolution (-gamma:2pi)
%
% NOTE:
% theta describes the trailing vane position (pala a seguire)
% following vane = "1", preceding vane = "2".
% vane is considered with its real geometry.

    %% DECLARATIONS AND PREALLOCATION %%
    R   = D/2;                      % stator radius [m]
    r   = d/2;                      % rotor radius [m]
    e   = R-r;                      % eccentricity [m]
    r_v = abs(r*sin(th_tilt));      % OA segment  [m]
    AB  = sqrt(r^2-r_v^2);          % [m]
    eta = sgn*(pi/2- abs(th_tilt)); % AOB angle [rad]
    OBB = r^2/2 *Gamma;             % circular sector [m2]

    %% VOLUMES CALCULATION %%
    V_cell       = NaN(size(theta)); % pre-allocation of 2D cell volume
    V_cell_limit = NaN(size(theta)); % pre-allocation of 1D cell volume

    % I PHASE: theta is between [-gamma:0)
    %          V_cell  : volume between tangency and leading vane [m3]
    %          alpha   : angle that define the tangency point A [rad]
    %          beta    : angle between  r and R in vane-stator contact point [rad]
    %          GAMMA   : compartement angular aperture seen from rotor center [rad]
    %          Ft      : factor used for 2D volume computation [m2]
    alpha1 = theta(1:Npt_cell) - eta;
    alpha2 = alpha1 + Gamma;
    AC2    = sgn*e*sin(alpha2)+sqrt(R^2-(r_v+e*cos(alpha2)).^2);
    OC2    = sqrt(AC2.^2+ r_v^2);
    BC2    = AC2-AB;
    phi2   = asin(AC2./OC2);
    psi2   = phi2-abs(eta);
    beta2  = asin((e*sin(pi-theta(1:Npt_cell)-Gamma-sgn*psi2))/R);
    GAMMA  = Gamma+theta(1:Npt_cell)+sgn*psi2-beta2;
    Ft1    = R^2/2*GAMMA;
    Ft2    = e*R/2*sin(GAMMA);
    Ft3    = sgn*(r_v/2*BC2);
    Ft4    = r^2/2*(Gamma+theta(1:Npt_cell));
    V_cell_limit(1:Npt_cell) = L*(Ft1 - Ft2 - Ft3 - Ft4);
    V_cell(1:Npt_cell) = V_cell_limit(1:Npt_cell) - Vol_Excluded(2,1:Npt_cell) - sgn*Vol_Excluded(3,1:Npt_cell);

    % II PHASE: theta is between [0:2pi-gamma)
    %           V_cell : volume between leading vane and trailing vane [m3]
    %           OCC    : circular sector OCC [m3]
    alpha1 = theta(pos_tang:end-Npt_cell) - eta;
    alpha2 = alpha1 +Gamma ;
    AC1    = sgn*e*sin(alpha1)+sqrt(R^2-(r_v+e*cos(alpha1)).^2);
    AC2    = sgn*e*sin(alpha2)+sqrt(R^2-(r_v+e*cos(alpha2)).^2);
    OC1    = sqrt(AC1.^2+ r_v^2);
    OC2    = sqrt(AC2.^2+ r_v^2);
    BC1    = AC1-AB;
    BC2    = AC2-AB;
    phi1   = asin(AC1./OC1);
    phi2   = asin(AC2./OC2);
    psi1   = phi1-abs(eta);
    psi2   = phi2-abs(eta);
    beta1  = asin((e*sin(pi-theta(pos_tang:end-Npt_cell)-sgn*psi1))/R);
    beta2  = asin((e*sin(pi-theta(pos_tang:end-Npt_cell)-Gamma-sgn*psi2))/R);
    THETA  = theta(pos_tang:end-Npt_cell)+sgn*psi1 -beta1;
    GAMMA  = beta1-beta2+Gamma-sgn*psi1+sgn*psi2;
    OCC    = R^2/2*GAMMA +e*R/2*sin(THETA)- e*R/2*sin(THETA+GAMMA);
    Ft1    = OCC - OBB;
    Ft2    = sgn*(r_v/2*BC1);
    Ft3    = sgn*(r_v/2*BC2);
    V_cell_limit(pos_tang:end-Npt_cell) = L*(Ft1 + Ft2 - Ft3);
    V_cell(pos_tang:end-Npt_cell) = V_cell_limit(pos_tang:end-Npt_cell) - Vol_Excluded(1,1:end-Npt_cell) + sgn*Vol_Excluded(3,1:end-Npt_cell) - Vol_Excluded(2,1+Npt_cell:end) - sgn*Vol_Excluded(3,1+Npt_cell:end);

    % III PHASE: theta is between [(2pi-gamma): 2pi]
    %            V_cell  : volume between trailing vane and tangency
    alpha1 = theta(end-Npt_cell+1:end) - eta;
    AC1    = sgn*e*sin(alpha1)+sqrt(R^2-(r_v+e*cos(alpha1)).^2);
    OC1    = sqrt(AC1.^2+ r_v^2);
    BC1    = AC1-AB;
    phi1   = asin(AC1./OC1);
    psi1   = phi1-abs(eta);
    beta1  = asin((e*sin(theta(end-Npt_cell+1:end)-pi+sgn*psi1))/R);
    GAMMA  = 2*pi-theta(end-Npt_cell+1:end)-sgn*psi1-beta1;
    Ft1    = R^2/2*GAMMA;
    Ft2    = e*R/2*sin(GAMMA);
    Ft3    = sgn*(r_v/2*BC1);
    Ft4    = r^2/2*(2*pi-theta(end-Npt_cell+1:end));
    V_cell_limit(end-Npt_cell+1:end) = L*(Ft1 - Ft2 + Ft3 - Ft4);
    V_cell(end-Npt_cell+1:end) = V_cell_limit(end-Npt_cell+1:end) - Vol_Excluded(1,end-Npt_cell+1:end) + sgn*Vol_Excluded(3,end-Npt_cell+1:end);

    % the first and last elements are set to zero
    V_cell(1)   = 0;
    V_cell(end) = 0;
    clear  alpha1 alpha2 AC1 AC2 OC1 OC2 BC1 BC2 phi1 phi2 beta1 beta2 psi1 psi2 GAMMA THETA OCC AB eta OBB Ft1 Ft2 Ft3 Ft4

    %% CHECKS %%
    % cell volume must be positive
    chck_1 = V_cell<0;
    if sum(chck_1) >= 1
        warning('S2_Volume_C:NegativeVolume','Negative cell volume for some angular position!');
        SX_Logfile ('w',{lastwarn});
    end

    % check on result in case of radial vanes
    if fDBG == 1
        if th_tilt == 0
            VcellOld        = NaN(size(theta));  % pre-allocation of 2D cell volume
            V_cellOld_limit = NaN(size(theta));  % pre-allocation of 1D cell volume

            % I PHASE: theta is between [-gamma:0)
            %         V_cellOLD  : volume between tangency and leading vane [m3]
            %         beta       : angle between  r and R in vane-stator contact point [rad]
            %         GAMMA      : compartement angular aperture seen from rotor center [rad]
            %         Ft         : factor used for 2D volume computation [m2]
            beta      = asin((e*sin(pi-theta(1:Npt_cell)-Gamma))/R);
            GAMMA     = theta(1:Npt_cell)+Gamma-beta;
            Ft1       = pi*(R^2-r^2);
            Ft2       = 0.5*R*(R*(2*pi-GAMMA)+e*sin(GAMMA));
            Ft3       = (0.5*r^2)*(2*pi-(theta(1:Npt_cell)+Gamma));
            V_cellOld_limit(1:Npt_cell) = L*(Ft1 - (Ft2 - Ft3));
            VcellOld(1:Npt_cell) = V_cellOld_limit(1:Npt_cell) - Vol_Excluded(2,1:Npt_cell) - sgn*Vol_Excluded(3,1:Npt_cell);

            % II PHASE: theta is between [0:2pi-gamma)
            %          V_cellOLD   : volume between leading vane and trailing vane [m3]
            %          alpha       : angle between r e R in the contact point for the trailing vane
            %          THETA       : trailing vane angular position as seen by stator
            alpha    = asin(e*sin(pi-theta(pos_tang:end-Npt_cell))/R);
            beta     = asin(e*sin(pi-theta(pos_tang:end-Npt_cell)-Gamma)/R);  % note pos_tang = N_vano+1
            THETA    = theta(pos_tang:end-Npt_cell)-alpha;
            GAMMA    = Gamma+theta(pos_tang:end-Npt_cell)-beta-THETA;
            Ft1      = pi*(R^2-r^2);
            Ft2      = 0.5*R*(THETA*R-e*sin(THETA));
            Ft3      = 0.5*r^2*theta(pos_tang:end-Npt_cell);
            Ft4      = 0.5*R*(R*(2*pi-THETA-GAMMA)+e*sin(THETA+GAMMA));
            Ft5      = (0.5*r^2)*(2*pi-theta(pos_tang:end-Npt_cell)-Gamma);
            V_cellOld_limit(pos_tang:end-Npt_cell) = L*(Ft1-(Ft2-Ft3)-(Ft4-Ft5));
            VcellOld(pos_tang:end-Npt_cell) = V_cellOld_limit(pos_tang:end-Npt_cell) - Vol_Excluded(1,1:end-Npt_cell) + sgn*Vol_Excluded(3,1:end-Npt_cell) - Vol_Excluded(2,1+Npt_cell:end) - sgn*Vol_Excluded(3,1+Npt_cell:end);

            % III PHASE: theta is between [(2pi-gamma): 2pi]
            %           V_cellOLD   : volume between trailing vane and tangency
            %           gammav      : cell angular aperture (vary from gamma to 0)
            gammav   = abs(theta(end-Npt_cell+1:end)-2*pi);
            beta     = asin(e*sin(pi-gammav)/R);
            GAMMA    = gammav-beta;
            Ft1      = pi*(R^2-r^2);
            Ft2      = 0.5*R*(R*(2*pi-GAMMA)+e*sin(GAMMA));
            Ft3      = (0.5*r^2)*(2*pi-gammav);
            V_cellOld_limit(end-Npt_cell+1:end) = L*(Ft1-(Ft2-Ft3));
            VcellOld(end-Npt_cell+1:end) = V_cellOld_limit(end-Npt_cell+1:end) - Vol_Excluded(1,end-Npt_cell+1:end) + sgn*Vol_Excluded(3,end-Npt_cell+1:end);

            % check on consistency of vanes
            err_mod = abs(V_cell - VcellOld);
            chck_2 = err_mod > toll_d;
            if sum(chck_2) >= 1
                warning('S2_Volume_C:Volume_cell','Volume of 1D and 2D models does not match')
                SX_Logfile ('d',{lastwarn});
            end

    if s == 0
        chck_3 = V_cell == V_cell_limit;
        chck_3 = sum(chck_3)/length(chck_3);
        if chck_3 ~= 1
            warning('Wrong volumes for 1D geometry')
            SX_Logfile ('w',{lastwarn});
        end
    end

            clear err_mod chck_1 chck_2 chck_3 s
        end
    end
end

