function [Fleak_g,Fleak_o,fOK] = S3_LeakageRS (MUmodel,RSmodel,c,process,D,d,L,RSclr,Gamma,TgAngle,dt_vano,theta,pos_SucOpen,pos_SucClose,pos_DisOpen,pos_DisClose,pos_End,V_cell,rpm,Pup,Pdo,T_mix,R_g,rho_g_RS,mu_g_RS,n_van,V_inj_nzl,cnd_suc,c_v,rho_o_RS,mu_o_RS,Npt_cell)
% This function compute leakages through clearance between stator and rotor sealing arc (path 1)
% 
% INPUT
% MUmodel [string]  : viscosity model for two-phase gas-liquid mixture
% RSmodel [string]  : rotor-stator leakage model
% c [-]             : stator geometry     1: circular   2: elliptical stator;
% process           : process selected (1: compression  2: expansion
% D [m]             : stator diameter
% d [m]             : rotor diameter
% L [m]             : rotor length
% RSclr [m]         : sealing arc clearance size
% Gamma [rad]       : angular extention of a cell (cavity between 2 vanes)
% TgAngle [rad]     : sealing arc angular extension
% dt_vano [s]       : time needed for a vane to sweep a cell [s]
% theta [rad]       : array of discretized angular position of trailing vane [-gamma:2pi]
% th_out [rad]      : updated delivery open angle (takes into account expansion ports flipping)
% pos_SucOpen [-]   : discretized suction open angle position in array theta
% pos_SucClose [-]  : discretized suction close angle position in array theta
% pos_DisOpen [-]   : discretized delivery open angle position in array theta
% pos_DisClose [-]  : discretized delivery close angle position in array theta
% pos_End [-]       : discretized position of cell end in array theta (is where tangency starts)
% V_cell [m^3]      : cell volumes array of a complete vane revolution (-gamma:2pi)
% rpm [rpm]         : rotor angular velocity
% Pup [pa]          : upstream pressure
% Pdo [pa]          : downstream pressure
% T_g [K]           : temperature array
% T_mix [K]         : oil + gas adiabatic mixing temperature for each discretization step
% R_g [J/kgK]       : specific gas constant
% rho_g_RS [kg/m3]  : leaked gas density
% mu_g_RS [Pa s]    : leaked gas dynamic viscosity
% n_van [-]         : number vanes
% V_inj_nzl [m3]    : volume of oil injected for each nozzle
% cnd_suc           : nozzles present during suction phase
% c_v [J/kg K]      : specific heat at costant volume
% rho_o_RS [kg/m3]  : leaked oil density
% mu_o_RS [Pa s]    : leaked oil kinematic viscosity
% Npt_cell [-]      : number of discretization points between 2 vanes

% OUTPUT
% Fleak_g [kg/s] : leaking gas mass flow rate through path 1
% Fleak_o [kg/s] : leaking oil mass flow rate through path 1

    %% DEFINITIONS %%
    % geometry
    r         = d/2;                        % rotor radius  [m]
    R         = D/2;                        % stator radius [m]
    e         = R-r;                        % eccentricity  [m]
    alpha     = 2*pi/c-theta(pos_DisClose); % angle between suction close and  minimum clearance position [0°] [rad]
   % beta     = alpha-TgAngle/2;            % angle between suction close and the beginning of tangency [rad]
    chi       = alpha+theta(pos_SucOpen);   % angle between delivery close and suction open [rad]
    rad_sec   = rpm/60*2*pi;                % velocity in rad per second [rad/s]
    Vtang_tot = RSclr*r*TgAngle*L;          % volume between rotor and stator tangency [m3]
    
    %% BULK PROPERTIES %%   
    if process == 1   % compressors
    V_g      = (V_cell(pos_DisOpen-Npt_cell)-sum(V_inj_nzl)); % gas volume before delivery opening [m3]
    Xvol_g   = V_g/V_cell(pos_DisOpen-Npt_cell);              % gas volumetric fraction before delivery open [m3 gas / (m3 gas + m3 oil)]
    mass_g   = V_g*rho_g_RS;                                  % gas mass before delivery opening [kg]
    mass_o   = sum(V_inj_nzl.*rho_o_RS);                      % oil mass before delivery opening [kg]
    mass_tot = sum([mass_g,mass_o],'omitnan');                % oil + gas mass before delivery opening [kg]
    Xmass_g  = mass_g/mass_tot;                               % gas mass fraction before delivery open [kg gas / (kg gas + kg oil)]
    
    else % expanders
    V_g      = (V_cell(pos_SucClose)-sum(V_inj_nzl(cnd_suc))); % gas volume before suction closing [m3]
    Xvol_g   = V_g/V_cell(pos_SucClose);                       % gas volumetric fraction before suction closing [m3 gas / (m3 gas + m3 oil)]    
    mass_g   = V_g*rho_g_RS;                                   % gas mass before suction closing [kg]
    mass_o   = sum(V_inj_nzl(cnd_suc).*rho_o_RS);              % oil mass before suction closing [kg]
    mass_tot = sum([mass_g,mass_o],'omitnan');                 % oil + gas mass before suction closing [kg]
    Xmass_g  = mass_g/mass_tot;                                % gas mass fraction before suction closing [kg gas / (kg gas + kg oil)]  
    end
    
    rho_b = 1/sum([Xmass_g/rho_g_RS,(1-Xmass_g)/rho_o_RS],'omitnan'); % bulk mixture density [kg/m3]
    T_b   = T_mix(end);                                               % bulk mixture temperature [kg/m3]
    
    switch MUmodel
        case 'Dukler'
            mu_b    = rho_b*sum([Xmass_g*mu_g_RS/rho_g_RS,(1-Xmass_g)*mu_o_RS/rho_o_RS],'omitnan');
        case 'Maxwell'
            mu_b_M1 = mu_o_RS.*(sum([2*mu_o_RS,mu_g_RS-2.*Xmass_g.*(mu_o_RS-mu_g_RS)],'omitnan')./sum([2*mu_o_RS,mu_g_RS+Xmass_g.*(mu_o_RS-mu_g_RS)],'omitnan'));            %-Maxwell 1 (gas based)
            mu_b_M2 = mu_g_RS.*(sum([2*mu_g_RS,mu_o_RS-2.*(1-Xmass_g).*(mu_g_RS-mu_o_RS)],'omitnan')./sum([2*mu_g_RS,mu_o_RS+(1-Xmass_g).*(mu_g_RS-mu_o_RS)],'omitnan'));    % -Maxwel 2 (liquid based)
            mu_b    = sum([mu_b_M1,mu_b_M2],'omitnan')/2; 
            if sum(isnan(mu_b_M1)) >= 1
                mu_b =  mu_b_M2;
            end
        case 'Awad'
            mu_b = 1/4*(sum([(3*Xmass_g-1)*mu_g_RS,(3*(1-Xmass_g)-1)*mu_o_RS],'omitnan') + sqrt(sum([sum([(3*Xmass_g-1)*mu_g_RS,(3*(1-Xmass_g)-1)*mu_o_RS],'omitnan')^2,8*mu_o_RS*mu_g_RS],'omitnan')));
        case 'MacAdams'
            mu_b = 1/sum([Xmass_g/mu_g_RS,(1-Xmass_g)/mu_o_RS],'omitnan');
        case 'Cicchitti'
            mu_b = sum ([Xmass_g*mu_g_RS,(1-Xmass_g)*mu_o_RS],'omitnan');
        case 'PureOil'
            mu_b = mu_o_RS;
    end
    
    %% MODEL SELECTION %%
    fOK  = 1;

    switch RSmodel
        %  Flow_Leak :  [kg/s] leaked mass flow rate when cell is communicating with port
        case 'Yanagisawa'
            FlowLeak = Yanagisawa (R,r,e,L,RSclr,TgAngle,Pup,Pdo,T_b,R_g,mu_b,n_van,c_v,Npt_cell);
        case 'Badr'
            FlowLeak = Badr (L,RSclr,Pup,Pdo,rho_b,T_b,R_g,c_v);
        case 'Poiselle-Couette'
            FlowLeak = PoiselleCouette (process,R,r,e,L,RSclr,rpm,TgAngle,Pup,Pdo,rho_b,mu_b,Gamma,Npt_cell);
        case 'Ishii'
            FlowLeak = Ishii (R,r,e,L,RSclr,TgAngle,Pup,Pdo,rho_b,mu_b,Gamma,Npt_cell);
        case 'Ferreira'
            FlowLeak = Ferreira (L,RSclr,Pup,Pdo,rho_b,mu_b);
        case 'Gasche'
            FlowLeak = Gasche (r,L,RSclr,Pup,Pdo,rho_b,mu_b,Xmass_g);
        otherwise
            fOK = 0;
            warning('S3_LeakagePath1:model','Invalid leakage model for rotor-stator clearance')
            SX_Logfile ('e',{lastwarn});
    end

    if process == 1 && Gamma >= alpha
        % The clearance can 'see' the delivery port. Leakage is divided in three phases
        % See Vallone 2018 for more info
        % PHASE I: trailing vane is between [-Gamma:Suction_Close]
        Dt_I     = (Gamma-alpha)/rad_sec;    % duration of Phase I [s]
        Mass_I_g = FlowLeak*Dt_I*Xvol_g;     % mass of gas leaking during phase I [kg]
        Mass_I_o = FlowLeak*Dt_I*(1-Xvol_g); % mass of oil leaking during phase I [kg]

        % PHASE II: trailing vane is between [Suction_Close:-TgAngle/2]
        %  !!! V_II must be corrected once the volume calculation will include
        %  the rotor slot in the stator (compenetrazione). Here I assumed that the volume loss due
        %  to the slot is equal to the cell volume at the beginning of the tangency
        %Dt_II     = beta/rad_sec;                         % duration of Phase II [s]
        V_II      = V_cell(pos_DisClose)- V_cell(pos_End); % cell volume at the beginning of phase II [m3]
        Mass_II_g = V_II*Xvol_g*rho_g_RS;                  % mass leaking during phase II [kg]
        Mass_II_o = V_II*(1-Xvol_g)*rho_o_RS;              % mass of oil leaking during phase II [kg]

        % PHASE III: trailing vane is between [-TgAngle:0]
        %Dt_III     = TgAngle/2/rad_sec;   % fraction of time when the vane is in tangency [s]
        Mass_III_o = Vtang_tot/2*rho_o_RS; % mass of leaking during phase III [kg]

        % TOTAL LEAKAGE COMPUTATION
        Mleak_g = (Mass_I_g + Mass_II_g);                             % mass of gas leaking during Dt_tot [kg]
        Mleak_o = sum([Mass_I_o , Mass_II_o , Mass_III_o],'omitnan'); % mass of oil leaking during Dt_tot [kg]
        Fleak_g = Mleak_g/dt_vano;                                     % leaking gas mass flow rate [kg/s]
        Fleak_o = Mleak_o/dt_vano;                                     % leaking oil mass flow rate [kg/s]
      
          
    elseif process == 1 && Gamma < alpha
        % the leakage mass is the one trapped in the cell after delivery closing. There is only phase II.
        % !!! V_Leak must be corrected once the volume alculation will include
        % the rotor slot in the stator.
        V_Leak   = V_cell(pos_DisClose);                                               % cell volume at the beginning of phase II [m3]
        Mleak_g  = V_Leak*Xvol_g*rho_g_RS;                                             % mass of gas leaking [kg]
        Mleak_o  = sum([V_Leak*(1-Xvol_g)*rho_o_RS , Vtang_tot/2*rho_o_RS],'omitnan'); % mass of oil leaking [kg]
        Fleak_g  = Mleak_g/(Gamma/rad_sec);                                            % leaking gas mass flow rate [kg/s]
        Fleak_o  = Mleak_o/(Gamma/rad_sec);                                            % leaking oil mass flow rate [kg/s]

            
    
    
        
        
        
    elseif process == 2
        % In expanders, due to mechanical reason, gas leakage is simply assumed to be equal to the one
        % given by models, only if suction and delivery ports are closer
        % than cell aperture angle. Otherwise no leakage occours.
        % NOTE: corrected by Franzetti-Persico according to Vallone (2017)
        % thesis.
        % =================================================================
        % =================================================================
%         if chi >= Gamma  % suction and delivery ports are too far
%             DTexp   = 0;
%         else             % ports close enough: leakage occours
%             DTexp   = (Gamma-chi)/rad_sec; % fraction of time when clearance can be the suction port
%         end

        DTexp   = (Gamma)/rad_sec; % fraction of time when clearance can be the suction port
        % =================================================================
        % =================================================================
        Fleak_g = FlowLeak*Xvol_g*DTexp/dt_vano;      % leaking gas mass flow rate [kg/s]
        Fleak_o = FlowLeak*(1-Xvol_g)*DTexp/dt_vano;  % mass of oil leaking [kg/s]
    end
    
end
                              
function [Mleak] = Yanagisawa (R,r,e,L,RSclr,TgAngle,Pup,Pdo,T_b,R_g,mu_b,n_van,c_v,Npt_cell)
%This function compute leakages through clearance between stator and rotor
%sealing arc (path 1), using Yanagisawa & Shimizu approach
% 
% INPUT
% R [m]          : stator diameter
% r [m]          : rotor diameter
% e [m]          : eccentricity
% L [m]          : rotor length
% RSclr [m]      : sealing arc clearance size
% TgAngle        : sealing arc angular extension
% Gamma [-]      : specific heat ratio
% Pup [pa]       : upstream pressure
% Pdo [pa]       : downstream pressure
% T_b [K]        : upstream bulk temperature   
% R_g [J/kgK]    : specific gas constant
% mu_b [Pa s]    : bulk dynamic viscosity
% n_van [-]      : number vanes
% % c_v [J/kg K] : specific heat at costant volume
% Npt_cell [-]   : number of discretization points between 2 vanes
% 
% OUTPUT
% Mleak [kg/s] : leakage mass flow rate through path 1
% 
% NOTES: - subscripts stand for:
%              up : upstream of clearance
%              do : downstream of clearance
%              t  : throat equivalernt channel
%              e  : equivalent channel exit
% HP: - Compresible flow
%     - Ideal gas
%     - Adiabatic Flow
%     - Viscous flow, either laminar or turbolent
%     - Fanno flow (entropy always increase. Max Mach is 1)
%
% Model inspired by: "Leakage losses with a rolling piston type compressor.
% I: radial clearanceon the rolling piston", Yanagisawa & Shimizu, 1985
%
% See Vallone 2018 for more info

%% DEFINITIONS
c_p   = c_v+R_g;  % gas specific heat at costant pressure [J/kgK]
Gamma = c_p/c_v;  % specific heat ratio [-]
options.Display = 'off';
options.MaxFunctionEvaluations = 1e5;
options.ConstraintTolerance = 1e-10;
options.FunctionTolerance = 1e-10;
options.StepTolerance     = 1e-10;

if TgAngle == 0 
    ArcLeak  = 2*pi/n_van/2;                                             % cavity angular extension [rad]
    phi      = linspace(-ArcLeak,ArcLeak,Npt_cell);                      % reference angle (equivalent to theta)
    R_th     = sqrt(R^2 -((e-RSclr).*sin(phi)).^2)-(e-RSclr).*cos(phi);  % stator radius seen from rotor center
    BCO      = asin(sin(phi)./R.*(e-RSclr));                             % stator center-stator surface-rotor center angle
    OBC      = pi-asin(sin(BCO)./r.*R_th);                               % rotor center-rotor surface-stator surface angle
    BOC      = pi-OBC-BCO;                                               % rotor surface-rotor center-stator surface angle
    Hc       = r./sin(BCO).*sin(BOC);                                    % equivalent channel height
    posRSclr = isnan(Hc);                                                % NaN values are spotted (they occour for phi=0)
    Hc(posRSclr) = RSclr;                                                % NaN values are changed
    Lf_0    = RSclr*trapz(R.*phi,1./Hc);  % equivalent channel specific length
    
else
    Lf_0 = r*TgAngle;
end

%% FIRST ATTEMPT SOLUTION: CHOKE HYPOTHESYS %%
% Mach number at exit is set to 1. Then it is verified if the hypothesis is confirmed
x_guess = [0.3, 0.1]';   % initial guess
x       = fsolve(@(x) Yanagawa_choke(x,T_b,L,RSclr,R_g,mu_b,Gamma,Lf_0,Pup),x_guess,options); % solve system

% solution estraction
Mt    = x(1);
Mleak = x(2);

% solve pressure ratios
Pt_e   = (1/Mt)*((Gamma+1)/(2+(Gamma-1)*Mt.^2))^0.5;
Pup_t  = (1+(Gamma-1)*Mt^2/2)^(Gamma/(Gamma-1));
Pratio = Pup_t*Pt_e;

%% COMPLETE SOLUTION %%
% if condition below is verified, flow is NOT choking, and the solution
% above is not valid. A new solution must be calculated.
if Pratio > (Pup/Pdo)   
    x_guess = [Mt, 0.8, 0.1]'; % initial guess
    x  = fsolve(@(x) Yanagawa_NOchoke(x,T_b,L,RSclr,R_g,mu_b,Gamma,Lf_0,Pup,Pdo),x_guess,options); % solve system
    
    % solution estraction
    Mleak = x(3);
end

end

function [Mleak] = Badr (L,RSclr,Pup,Pdo,rho,T_b,R_g,c_v)
%This function compute leakages through clearance between stator and rotor
%sealing arc (path 1), using Badr compressible approach
% 
% INPUT
% L [m]         : rotor length
% RSclr [m]     : sealing arc clearance size
% Pup [pa]      : upstream pressure
% Pdo [pa]      : downstream pressure
% rho [kg/m3]   : upstream gas density
% T_b  [K]      : upstream bulk temperature
% R_g [J/kgK]   : specific gas constant
% c_v [J/kgK]   : specific heat at costant volume
% 
% OUTPUT
% Mleak [kg/s] : leakage mass flow rate through path 1
%  
% HP: - Comnpressible flow
%     - inviscid flow
%     - isentropic nozzle
%
% Model inspired by: "Multi-Vane Expanders: Internal-Leakage Losses", O. Badr, 1985
% See Vallone 2018 thesis for more info

%% DEFINITIONS %%
Cd    = 1;        % discharge coefficient
c_p   = c_v+R_g;  % gas specific heat at costant pressure [J/kgK]
Gamma = c_p/c_v;  % specific heat ratio [-]

if Pdo/Pup < (2/(Gamma+1))^(Gamma/(Gamma-1))
    Pt = Pup*(2/(Gamma+1))^(Gamma/(Gamma-1));
else
    Pt = Pdo;
end

Mleak = Cd*RSclr*L*rho*sqrt(2*c_p*T_b*((Pt/Pup)^(2/Gamma)-(Pt/Pup)^((Gamma+1)/Gamma)));
end

function [Mleak] = PoiselleCouette (process,R,r,e,L,RSclr,rpm,TgAngle,Pup,Pdo,rho_b,mu_b,Gamma,Npt_cell)
%This function compute leakages through clearance between stator and rotor
%sealing arc (path 1), using poiselle-Couette approach for 2D flow
% 
% INPUT
% process      : process selected (1: compression  2: expansion
% R [m]        : stator diameter
% r [m]        : rotor diameter
% e [m]        : eccentricity
% L [m]        : rotor length
% RSclr [m]    : sealing arc clearance size
% rpm [rpm]    : rotor velocity
% TgAngle      : sealing arc angular extension
% Pup [pa]     : upstream pressure
% Pdo [pa]     : downstream pressure
% rho_b        : upstream bulk density
% mu_b [Pa s]  : bulk dynamic viscosity
% Gamma [rad]  : angular extention of a cell (cavity between 2 vanes)
% Npt_cell [-] : number of discretization points between 2 vanes
% 
% OUTPUT
% Mleak [kg/s] : leakage mass flow rate through path 1
%  
% HP: - Uncompresible flow
%     - Laminar flow
%     - Viscous flow
%     - Poiselle-Covette flow (wall velocity accounted)
%
% Model inspired by: "Simulation and comparison of leakage characteristics of R290 in rolling piston type rotary compressor", Dehua Cai, 2015
% See Vallone 2018 thesis for more info

    %% DEFINITIONS %%
    Uwall = rpm/60*2*pi; % rotor surface velocity [rad/s]

    if process == 1 % compressione
        s = 1;
    else            % espansione
        s = -1;
    end
    
    %% CALCULATIONS %%
    if TgAngle == 0
        phi     = linspace(-Gamma,Gamma,Npt_cell);                          % angular extention of leakage model [rad]
        R_th    = sqrt(R^2 -((e-RSclr).*sin(phi)).^2)-(e-RSclr).*cos(phi);  % stator radius from rotor center [m]
        delta   = R_th - r;                                                 % clearance size, function of phi [m]
        Leq     = trapz(phi.*r,(1./delta).^3);                              % fictional sealing arc extension [rad]
    else
        Leq = TgAngle*r/RSclr^3;                                            % real sealing arc extension [rad]
    end

        Mleak_P = 1/(12*mu_b*Leq)*L*(Pup-Pdo)*rho_b;  % leakage mass flow rate due to pressure gradient [kg/s]
        Mleak_V = L*rho_b*RSclr/2*r*Uwall;            % leakage mass flow rate due to wall velocity [kg/s]
        Mleak   = Mleak_P + s * Mleak_V;              % total leakage mass flow rate through path 1 [kg/s]
end

function [Mleak] = Ishii (R,r,e,L,RSclr,TgAngle,Pup,Pdo,rho_b,mu_b,Gamma,Npt_cell)
%This function compute leakages through clearance between stator and rotor
%sealing arc (path 1), using Ishii aproach.
% 
% INPUT
% R [m]        : stator diameter
% r [m]        : rotor diameter
% e [m]        : eccentricity
% L [m]        : rotor length
% RSclr [m]    : sealing arc clearance size
% TgAngle      : sealing arc angular extension
% Pup [pa]     : upstream pressure
% Pdo [pa]     : downstream pressure
% mu_b [Pa s]  : dynamic viscosity
% rho_b        : upstream density
% Gamma [rad]  : angular extention of a cell (cavity between 2 vanes)
% Npt_cell [-] : number of discretization points between 2 vanes
% 
% OUTPUT
% Mleak [kg/s] : leakage mass flow rate through path 1
%  
% HP: - Uncompresible flow
%     - Turbolent flow
%     - Viscous flow
%
% Model inspired by: "Refrigerant Leakage Flow Evaluation for Scroll Compressors”", N. Ishii, 1996
% See Vallone 2018 thesis for more info
    
    %% CALCULATIONS %%
    options.Display = 'off';
    if TgAngle == 0
        phi     = linspace(-Gamma,Gamma,Npt_cell);                          % angular extention of leakage model [rad]
        R_th    = sqrt(R^2 -((e-RSclr).*sin(phi)).^2)-(e-RSclr).*cos(phi);  % stator radius from rotor center [m]
        delta   = R_th - r;                                                 % clearance size, function of phi [m]
        Deq     = delta.*L./(2.*(delta + L));                               % hydraulic diameter of clearance
        Leq     = phi.*r;                                                   % fictional sealing length
        FUN     = @(Uo) trapz(Leq,0.35.*(4.*Deq.*Uo.*RSclr./delta.*rho_b./mu_b).^(-0.25) .*(Uo.*RSclr./delta).^2 ./ (8.*Deq)) - (Pup-Pdo)./rho_b;
        Uo      = fsolve(FUN, 0.1, options);
    else
        Deq = RSclr.*L/(2.*(RSclr + L));  % real sealing arc extension [rad]
        Leq = TgAngle.*r;                 % sealing length [m]
        FUN = @(Uo) 0.35*(4*Deq*Uo.*rho_b/mu_b)^(-0.25)*Leq.*Uo^2/(8.*Deq) - (Pup-Pdo)/rho_b; 
        Uo  = fsolve(FUN, 0.1,options);
    end
    
    Mleak = Uo*RSclr*L*rho_b; % leakage mass flow rate through path 1 [kg/s]
end
    
function [Mleak] = Ferreira (L,RSclr,Pup,Pdo,rho_b,mu_b)
%This function compute leakages through clearance between stator and rotor
%sealing arc (path 1), using Ferreira uncompressible approach
% 
% INPUT
% L [m]         : rotor length
% RSclr [m]     : sealing arc clearance size
% Pup [pa]      : upstream pressure
% Pdo [pa]      : downstream pressure
% rho_b [kg/m3] : upstream bulk density
% mu_b [Pa s]   : bulk dynamic viscosity
% 
% OUTPUT
% Mleak [kg/s] : leakage mass flow rate through path 1
%  
% HP: - semiempirical correlation
%     - based on poiselle equation
%
% Model inspired by: "Bicylindrical Coordinate Formulation for the Leakage Flow Through the Minimal Clearance in a Rolling Piston Compressor", R.T. Ferreira, 1992
% See Vallone 2018 thesis for more info
%% DEFINITIONS %%
Mleak = 0.0162*mu_b*RSclr*(Pup-Pdo)*(rho_b*RSclr^2/mu_b^2)*(L/RSclr)^0.504;
end

function [Mleak] = Gasche (r,L,RSclr,Pup,Pdo,rho_b,mu_b,Xmass_g)
%This function compute leakages through clearance between stator and rotor
%sealing arc (path 1), using Gasche uncompressible mix approach
% 
% INPUT
% r [m]         : rotor diameter
% L [m]         : rotor length
% RSclr [m]     : sealing arc clearance size
% Pup [pa]      : upstream pressure
% Pdo [pa]      : downstream pressure
% rho_b [kg/m3] : upstream density
% mu_b [Pa s]   : dynamic viscosity
% Xmass_g [-]   : gas mass fraction [kg gas/kg mix]
% 
% OUTPUT
% Mleak [kg/s] : leakage mass flow rate through path 1
%  
% HP: - semiempirical correlation
%     - based on poiselle equation
%
% Model inspired by: "A model to predict R134a refrigerant leakage through the radial clearance of a rolling piston comprssor", J.L.. Gasche, 2012
% See Vallone 2018 thesis for more info
%% DEFINITIONS %%
XmassLim = 0.2211/0.276;
if Xmass_g>XmassLim
    Xmass_g = XmassLim;
    warning('S3_LeakageRS')
    warning('S3_LeakageRS:Gasche','LEAK.RSmodel:Gas mass fraction exceed maximum for Gasche model');
    SX_Logfile ('w',{lastwarn});
end
Mleak = 1/12*(0.2211-0.276*Xmass_g)*rho_b*L*(Pup-Pdo)*RSclr^2.5/mu_b/r^0.5;
end

function FUN = Yanagawa_choke(x,Tup,L,RSclr,R_g,mu_b,Gamma,Lf_0,Pup)
% This function set the system of 2 non linear equation under choking hypothesys
% Unknowns of the system:
% Mt = x(1)            : throat mach number [-] 
% Mleak = x(2)     : leakage mass flow rate [kg/s]

%% PREAMBLE %%
FUN = NaN(2,1);
Me  = 1;        

% define lambda
Re     = 2*x(2)/(mu_b*L);
lambda = 96/Re*(Re<=3560) + 0.3164/Re^0.25*(Re>3560);

%% SYSTEM DEFINITION %%
FUN (1) = abs((1-x(1)^2)/(Gamma*x(1)^2)+(Gamma+1)/(2*Gamma)*log((Gamma+1)*x(1)^2/(2+(Gamma-1)*x(1)^2))/(Lf_0*lambda/(2*RSclr)))-1; % I equation

% Pressure ratios
Pt_e   = (1/x(1))*((Gamma+1)/(2+(Gamma-1)*x(1).^2))^0.5; % ratio between throat pressure and exit pressure [-]
Pup_t  = (1+(Gamma-1)*x(1)^2/2)^(Gamma/(Gamma-1));       % ratio between upstream pressure and throat pressure [-]
Pratio = Pup_t*Pt_e;                                     % ratio between upstream and exit pressure [-]

% Define exit condition
Te = Tup/(1+Me^2*(Gamma-1)/2);  % temperature at channel exit [K]
Pe = Pup/Pratio;                % pressure  at channel exit [Pa]
Ve = Me*(Gamma*R_g*Te)^0.5;     % velocity at channel exit [m/s]

% Mass flow rate 
Mleak  = RSclr*L*Ve*Pe/(R_g*Te);  % leakage mass flow rate [kg/s]
FUN(2) = Mleak-x(2);              % II equation
end


function FUN = Yanagawa_NOchoke(x,Tup,L,RSclr,R_g,mu_b,Gamma,Lf_0,Pup,Pdo)
% This function set the system of 3 non linear equation for not-choking flow
% Unknowns of the system:
% Mt = x(1)         : throat mach number [-]
% Me = x(2)         : equivalent channel exit mach number [-]
% Mleak = x(3)  : leakage mass flow rate [kg/s]

%% PREAMBLE %%
FUN = NaN(3,1);

% define lambda
Re     = 2*x(3)/(mu_b*L);
lambda = 96/Re*(Re<=3560) + 0.3164/Re^0.25*(Re>3560);

%% SYSTEM DEFINITION %%
Lf_i   = (2*RSclr)/lambda*((1-x(1)^2)/(Gamma*x(1)^2)+(Gamma+1)/(2*Gamma)*log((Gamma+1)*x(1)^2/(2+(Gamma-1)*x(1)^2))); % new specific length [m]
FUN(1) =                   (1-x(2)^2)/(Gamma*x(2)^2)+(Gamma+1)/(2*Gamma)*log((Gamma+1)*x(2)^2/(2+(Gamma-1)*x(2)^2))-(Lf_i-Lf_0)*lambda/(2*RSclr);  % I equation

% pressure ratios
Pt_x   = (1/x(1))*((Gamma+1)/(2+(Gamma-1)*x(1).^2))^0.5; % ratio between throat pressure and local pressure [-]
Pup_t  = (1+(Gamma-1)/2*x(1)^2)^(Gamma/(Gamma-1));       % ratio between upstream pressure and throat pressure [-]
Pe_x   = (1/x(2))*((Gamma+1)/(2+(Gamma-1)*x(2).^2))^0.5; % ratio between exit pressure and local pressure [-]
Pratio = Pup_t*Pt_x/Pe_x;                                % ratio between upstream and exit pressure [-]
FUN(2) = abs(Pratio/(Pup/Pdo))-1;  % II equation

% Define exit condition
Te = Tup/(1+x(2)^2*(Gamma-1)/2); % temperature at channel exit [K]
Pe = Pup/Pratio;                 % pressure  at channel exit [Pa]
Ve = x(2)*(Gamma*R_g*Te)^0.5;    % velocity at channel exit [m/s]

% Mass flow rate 
Mleak  = RSclr*L*Ve*Pe/(R_g*Te);  % leakage mass flow rate [kg/s]
FUN(3) = abs(Mleak /x(3))-1;      % III equation
end

