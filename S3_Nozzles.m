function [V_nz_tot,m_nz_tot,SMD,D_g,f_g,fOK] = S3_Nozzles(NOZZLES,n,rho_x,nu_x,mu_x,sigma_x,p_g,T_g,R_g)
% This function characterizes the nozzles evaluating volume and mass flow rate of oil,
% Sauter Mean Diameter, mean droplet diameter for each class and class frequency.
%
% INPUT
% NOZLES             : NOZZLE structure
% rho_x [kg/m3]      : oil density at injection temperature
% nu_x [cSt]         : oil kinematic viscosity at injection temperature
% mu_x [Pa*s]        : oil dynamic viscosity at injection temperature
% sigma_x [N/m]      : oil surface tension at injection temperature
% p_inj [Pa]         : oil injection pressure
% p_g [Pa]           : chamber gas pressure
% T_g [K]            : chamber gas temperature
% R_g [J/kgK]        : specific gas constant
%
% OUPUT
% V_nz_tot [m3/s] : total oil volume flow rate through the nozzles
% m_nz_tot [kg/s] : total oil mass flow rate through the nozzles
% SMD [m]         : Sauter Mean Diameter
% D_g [m]         : mean droplet diameter for each class
% f_g [adim]      : frequency for each class
% fOK [adim]      : OK flag (=0 if a problem occurs, =1 if no problems occur)

    %% DECLARATIONS %%
    fOK = 1; % OK flag (=0 if a problem occurs, =1 if no problems occur)

    %% UNPACKING %%
    type_nz  = NOZZLES.type_nz{n};  % nozzle type [string]
    num_nz   = NOZZLES.num_nz(n);   % number of nozzle [-]
    num_cls  = NOZZLES.num_cls(n);  % number of classes [-]
    alpha    = NOZZLES.alpha_nz(n); % probability that droplet diameters fall below d_min (tail at alpha/2 of probability) and above d_max (tail at 1-alpha/2)
    p_inj    = NOZZLES.p_inj(n);    % injection pressurs [Pa]
    
    %% MAIN CALCULATIONS %%
    % known flow nozzles
    if strcmp(type_nz,'kf')
        D_g      = NOZZLES.specs{n}.D_g;           % class droplet diameter, from user input [m]
        f_g      = NOZZLES.specs{n}.f_g;           % class frequency [-]
        V_nz_tot = NOZZLES.specs{n}.V_inj*num_nz;  % total volume flow rate through the nozzles [m3/s]
        m_nz_tot = V_nz_tot*rho_x(n);              % total oil mass flow rate through the nozzles [kg/s]
        SMD      = sum(D_g.*f_g);

    % pressure driven nozzle (hc, fc)
    else
        Delta_p = p_inj-p_g;                  % pressure difference across the injector [Pa]
        rho_g   = p_g/(R_g*T_g);              % ideal gas density [kg/m3] !!! TO BE MODIFIED IN FUTURE !!!
        D_or    = NOZZLES.specs{n}.D_or/1000; % orefice diameter[m]
        A_or    = pi*(D_or)^2/4;              % discharge orifice area [m2]
        u       = sqrt(2*(Delta_p)/rho_x);    % ideal mean injection velocity [m/s]
        Re      = (rho_x*u*D_or)/mu_x;        % Reynolds Number
        We      = (rho_g*u^2*D_or)/sigma_x;   % Weber Number

        if Delta_p <= 0
            % if injection pressure is below cell pressure, no oil can be
            % injected; All output variables are set to NaN and fOK is set to zero
            fOK      = 0;
            V_nz_tot = NaN;
            m_nz_tot = NaN;
            SMD      = NaN;
            D_g      = NaN(1,num_cls);
            f_g      = NaN(1,num_cls);
            warning('S3_Nozzles:pressure','Nozzles ''%s'' injection pressure is lower than chamber pressure', type_nz);
            SX_Logfile ('e',{lastwarn});
        end

        % plain orifice nozzles
        if strncmp(type_nz,'po',2) && fOK
            L_or = NOZZLES.specs{n}.L_or/1000; % orifice lenght
            [V_nz_tot,m_nz_tot,Cd] = flowrate_nozzle_plain_orifice(num_nz,rho_x,D_or,L_or,A_or,u,Re); % volumetric and mass flow rates [L/min, kg/s]
            switch type_nz
                case 'po'
                    [SMD] = SMD_nozzle_plain_orifice(nu_x,D_or,u);                         %Merrington e Richardson
                case 'po2'
                    [SMD] = SMD_nozzle_plain_orifice2(rho_x,mu_x,sigma_x,rho_g,D_or,u,Cd); %Varde
                case 'po3'
                    [SMD] = SMD_nozzle_plain_orifice3(rho_x,nu_x,sigma_x,Delta_p,rho_g);   %Elkotb
                case 'po4'
                    [SMD] = SMD_nozzle_plain_orifice4(rho_g,Delta_p,num_nz,V_nz_tot);      %Hiroyasu e Katoda
            end
            [D_g, f_g] = SizeDistribution(type_nz,num_cls,alpha,SMD/6);

        % hollow cone nozzles
        elseif strncmp(type_nz,'hc',2) && fOK
            gamma_nz            = NOZZLES.specs{n}.gamma_nz;    % angular aperture of cone [°]
            V_ref               = NOZZLES.specs{n}.V_ref;       % water flow rate supplied with the reference pressure drop [l/min]
            Delta_p_ref         = NOZZLES.specs{n}.Delta_p_ref; % reference pressure drop [Pa] 
            [V_nz_tot,m_nz_tot] = flowrate_nozzle_pressure_swirl(num_nz,rho_x,Delta_p,V_ref,Delta_p_ref,type_nz);         % volumetric and mass flow rates [L/min, kg/s]
            [SMD]               = SMD_nozzle_cone_hollow(m_nz_tot,num_nz,rho_x,mu_x,sigma_x,Delta_p,rho_g,D_or,gamma_nz); % sauter mean diameter [m]
            [D_g, f_g]          = SizeDistribution(type_nz,num_cls,alpha,SMD/6); % droplet diameter and frequency [m, -]

        % full cone nozzles
        elseif strncmp(type_nz,'fc',2) && fOK
            V_ref               = NOZZLES.specs{n}.V_ref;       % water flow rate supplied with the reference pressure drop [l/min]  
            Delta_p_ref         = NOZZLES.specs{n}.Delta_p_ref; % reference pressure drop [Pa] 
            [V_nz_tot,m_nz_tot] = flowrate_nozzle_pressure_swirl(num_nz,rho_x,Delta_p,V_ref,Delta_p_ref,type_nz); % volumetric and mass flow rates [L/min, kg/s]
            switch type_nz
                case 'fc'
                    [SMD] = SMD_nozzle_cone_full(D_or,Re,We); %Estes and Mudawar
                case 'fc2'
                    [SMD] = SMD_nozzle_cone_full2(rho_x,mu_x,sigma_x,Delta_p,rho_g,D_or); % Kohnen et al.
                case 'fc3'
                    [SMD] = SMD_nozzle_cone_full(D_or,Re,We); %Estes and Mudawar with SMD reduced by half
                    SMD  = SMD/2;
            end
            [D_g, f_g] = SizeDistribution(type_nz,num_cls,alpha,SMD/6); % droplet diameter and frequency [m, -]
        end
    end
end    
 
%% AUXILIARY FUNCTIONS %%
    
function [D_g, f_g] = SizeDistribution(type_nz,num_cls,alpha,SMD)
% INPUT 
% type_nz [string] : nozzle type
% num_cls [-]      : number of class
% alpha [-]        : probability that droplet diameters fall below d_min
%                        (tail at alpha/2 of probability) and above d_max (tail at 1-alpha/2)
% SMD [m]          : sauter mean diameter
% 
% OUTPUT
% D_g [m]    : mean droplet diameter for each class
% f_g [adim] : frequency for each class

% NOTE: Droplet size distribution
%       The Rosin Rammler volume distribution equation is adopted here for the droplet size.
%       Rosin Rammler probability density function     : f(d) = beta/delta*(d/delta)^(beta-1)*exp(-(d/delta)^beta)
%       Rosin Rammler cumulative distribution function : F(d) = 1- exp(-(d/delta)^beta)

    %% CALCULATIONS %%   
    beta  = 3;                     % distribution parameter
    delta = SMD*gamma(1-1/beta);   % size parameter
    mu_RR = delta*gamma(1+1/beta); % mean of the Rosin Rammler distribution

    if strncmp(type_nz,'po',2) % Plain Orifice
        D_g = mu_RR;
        f_g = 1;
    else  % Hollow cone and full cone
        if num_cls == 1
            D_g = mu_RR;
            f_g = 1;
        else
            [D_g,f_g] = discret_RR(num_cls,delta,beta,alpha);
        end
    end
end

function [V_nz_tot,m_nz_tot,Cd] = flowrate_nozzle_plain_orifice(num_nz,rho_x,D_or,L_or,A_or,u,Re)
%This function calculates volume and mass flow rate of oil for plain orifices.
%
%INPUT
%num_nz [adim]  : number of plain orifices
%rho_x [kg/m3]  : oil density at injection temperature
%D_or [m]       : orifice diameter
%L_or [m]       : orifice length
%A_or [m2]      : discharge orifice area
%u [m/s]        : ideal mean injection velocity
%Re [adim]      : Reynolds number
%
%OUTPUT
%V_nz_tot [m3/s] : total volume flow rate through the orifices
%m_nz_tot [kg/s]  : total mass flow rate through the orifices
%
%EXAMPLE
%[V_nz_tot,m_nz_tot,Cd] = flowrate_nozzle_plain_orifice(11,905,10,4,1.2568e-05,21.0239,1001.3836)

%Cdu and Cd are calculated according correlations proposed in
%"Discharge Coefficients for Incompressible Non-Cavitating Flow Through
%Long Orifices" (Lichtarowicz et al, 1965)
    Cdu  = 0.827-0.0085*L_or/D_or;             % ultimate discharge coefficient
    Cd   = 1/(1/Cdu+20*(1+2.25*L_or/D_or)/Re); % discharge coefficient
    m_nz = rho_x*u*A_or*Cd;                    % mass flow rate for single orifice [kg/s]
    V_nz = m_nz/rho_x;                         % volume flow rate for single orifice [m3/s]
    %all nozzles
    m_nz_tot = m_nz*num_nz;
    V_nz_tot = V_nz*num_nz;
end

function [SMD] = SMD_nozzle_plain_orifice(nu_x,D_or,u)
%This function calculates the Sauter Mean Diameter of plain orifices according the correlation
%proposed in "The Break-up of liquid jets" (Merrington and Richardson, 1947)
%
%INPUT
%nu_x [cSt] : oil kinematic viscosity at injection temperature
%D_or [m]   : orifice diameter
%u [m/s]    : mean injection velocity
%
%OUTPUT
%SMD [m]    : Sauter Mean Diameter
%
%EXAMPLE
%[SMD] = nozzle_plain_orifice3(63.3804,4,21.0239)

    nu_x = nu_x/10^6;               %[m^2/s]
    SMD  = 500*D_or^1.2*nu_x^0.2/u; %[m]
end

function [SMD] = SMD_nozzle_plain_orifice2(rho_x,mu_x,sigma_x,rho_g,D_or,u,Cd)
%This function calculates the Sauter Mean Diameter of plain orifices according the correlation
%proposed in "Spray Angle and Atomization in Diesel Spray" (Varde et al, 1984)
%
%INPUT
%rho_x [kg/m3] : oil density at injection temperature
%mu_x [cSt]    : oil dynamic viscosity at injection temperature
%sigma_x [N/m] : oil surface tension at injection temperature
%rho_g [kg/m3] : ideal gas density
%D_or [m]      : orifice diameter
%u [m/s]       : ideal mean injection velocity
%C_d [adim]    : discharge coefficient
%
%OUTPUT
%SMD [m]       : Sauter Mean Diameter
%
%EXAMPLE
%[SMD] = SMD_nozzle_plain_orifice2(905,0.0760,0.0323,2.1557,4,21.0239,0.65)

    u_real  = u*Cd; %real mean injection velocity [m/s]
    Re_real = (rho_x*u_real  *D_or)/mu_x;
    We_real = (rho_g*u_real^2*D_or)/sigma_x;
    SMD     = D_or*8.7*(Re_real*We_real)^(-0.28); %[m]
end

function [SMD] = SMD_nozzle_plain_orifice3(rho_x,nu_x,sigma_x,Delta_p,rho_g)
%This function calculates the Sauter Mean Diameter of plain orifices via the correlation
%proposed in "Fuel Atomization for Spray Modeling" (Elkotb, 1982)
%
%INPUT
%rho_x [kg/m3] : oil density at injection temperature
%nu_x [cSt]    : oil kinematic viscosity at injection temperature
%sigma_x [N/m] : oil surface tension at injection temperature
%Delta_p [Pa]  : pressure difference across the injector
%rho_g [kg/m3] : ideal gas density
%
%OUTPUT
%SMD [m]       : Sauter Mean Diameter
%
%EXAMPLE
%[SMD] = SMD_nozzle_plain_orifice(905,63.3804,0.0323,5*10^5,2.1557)

    nu_x    = nu_x/10^6;    %[m^2/s]
    Delta_p = Delta_p/10^5; %[bar]
    SMD     = 6156*nu_x^0.385*sigma_x^0.737*rho_x^0.737*rho_g^0.06*(Delta_p)^(-0.54); %[10^-6 m]
    SMD     = SMD/10^6;     %[m]
end

function [SMD] = SMD_nozzle_plain_orifice4(rho_g,Delta_p,num_nz,V_nz_tot)
%This function calculates the Sauter Mean Diameter of plain orifices according the correlation
%proposed in "Droplet Size Distribution in a Diesel Combustion Chamber"
%(Hiroyasu and Katoda, 1974).
%
%INPUT
%rho_g [kg/m3]   : ideal gas density
%Delta_p [Pa]    : pressure difference across the injector
%num_nz [adim]   : number of nozzles
%V_nz_tot [m3/s] : total volume flow rate through the orifices
%
%OUTPUT
%SMD [m]         : Sauter Mean Diameter
%
%EXAMPLE
%[SMD] = nozzle_plain_orifice4(2.1557,5*10^5,11,215.4127)

    Q   = V_nz_tot/num_nz; %volume flow rate for single nozzle [m3/s]
    SMD = 2330*rho_g^0.121*Q^0.131*Delta_p^(-0.135); %[10^-6*m]
    SMD = SMD/10^6; %[m]
end

function [V_nz_tot,m_nz_tot] = flowrate_nozzle_pressure_swirl(num_nz,rho_x,Delta_p,V_ref,Delta_p_ref,type_nz)
%This function calculates volume and mass flow rate of oil for pressure-swirl nozzles.
%
%INPUT
%num_nz [adim]    : number of nozzles
%rho_x [kg/m3]    : oil density at injection temperature
%Delta_p [Pa]     : pressure difference across the nozzle
%V_ref [l/min]    : reference volume flow rate for the nozzle
%Delta_p_ref [Pa] : reference pressure difference across the nozzle
%type_nz [string] : injector type
%
%OUTPUT
%V_nz_tot [m3/s] : total volume flow rate through the nozzles
%m_nz_tot [kg/s]  : total mass flow rate through the nozzles
%
%EXAMPLE
%[V_nz_tot,m_nz_tot] = flowrate_nozzle_pressure_swirl(2,905,5*10^5,12.5,15*10^5,0.4)

%Flow rate calculations
%Volume flow rate through the nozzle is derived from a correlation
%provided by the manufacturer
    F    = sqrt(1000/rho_x);                    % density conversion factor
    if strncmp(type_nz,'hc',2)
        V_nz = V_ref*(Delta_p/Delta_p_ref)^0.5; % water volume flow rate for single nozzle [l/min]
    elseif strncmp(type_nz,'fc',2)
        V_nz = V_ref*(Delta_p/Delta_p_ref)^0.4; % water volume flow rate for single nozzle [l/min]
    end
    V_nz = F*V_nz/1000/60; % water to oil conversion [m3/s]
    m_nz = V_nz*rho_x;     % mass flow rate for single nozzle [kg/s]                             
    % all nozzles
    V_nz_tot = V_nz*num_nz;
    m_nz_tot = m_nz*num_nz;
end

function [SMD] = SMD_nozzle_cone_hollow(m_nz_tot,num_nz,rho_x,mu_x,sigma_x,Delta_p,rho_g,D_or,gamma_nz)
%This function calculates Sauter Mean Diameter for hollow cone pressure-swirl nozzles
%via the correlation proposed in "Mean Drop Sizes from Pressure-Swirl Nozzles" (Wang and Lefebvre, 1987)
%
%INPUT
%m_nz_tot [kg/s] : total mass flow rate through the nozzles
%num_nz [adim]   : number of nozzles
%rho_x [kg/m3]   : oil density at injection temperature
%mu_x [Pa*s]     : oil dynamic viscosity at injection temperature
%sigma_x [N/m]   : oil surface tension at injection temperature
%Delta_p [Pa]    : pressure difference across the nozzle
%rho_g [kg/m3]   : gas density
%D_or [m]        : orifice diameter
%gamma_nz [deg]  : spray cone angle of the nozzle
%
%OUTPUT
%SMD [m]         : Sauter Mean Diameter
%
%EXAMPLE
%[SMD] = SMD_nozzle_cone_hollow(0.0448,2,905,0.0760,0.0323,5*10^5,2.1557,8,45)

    gamma_nz = deg2rad(gamma_nz); %[rad]
    m_nz     = m_nz_tot/num_nz;   %mass flow rate for single nozzle [kg/s]
    t_costh  = 2.7*((D_or*m_nz*mu_x)/(rho_x*Delta_p))^0.25*cos(gamma_nz/2); % lubricant film thickness in final orifice [m]
    SMD      = 4.52*((sigma_x*mu_x^2)/(rho_g*Delta_p^2))^0.25*(t_costh)^0.25...
               + 0.39*((sigma_x*rho_x )/(rho_g*Delta_p  ))^0.25*(t_costh)^0.75;
end

function [SMD] = SMD_nozzle_cone_full(D_or,Re,We)
%This function calculates Sauter Mean Diameter for full cone pressure-swirl nozzles
%via the correlation proposed in "Correlation of Sauter mean diameter and critical
%heat flux for spray cooling of small surfaces" (Estes and Mudawar, 1995)
%
%INPUT
%D_or [m]  : orifice diameter
%Re [adim] : Reynolds number
%We [adim] : Weber number
%
%OUTPUT
%SMD [m]   : Sauter Mean Diameter
%
%EXAMPLE
%[SMD] = SMD_nozzle_cone_full(10,1819.3553,204.7176)

    SMD = 3.67*D_or*(Re*We^0.5)^(-0.259);
end

function [SMD] = SMD_nozzle_cone_full2(rho_x,mu_x,sigma_x,Delta_p,rho_g,D_or)
%This function calculates Sauter Mean Diameter for full cone pressure-swirl nozzles
%via the correlation introduced in "Measurement of the droplet size of a full cone nozzle"
%(Kohnen et al., 2010)
%
%INPUT
%rho_x [kg/m3] : oil density at injection temperature
%mu_x [Pa*s]   : oil dynamic viscosity at injection temperature
%sigma_x [N/m] : oil surface tension at injection temperature
%Delta_p [Pa]  : pressure difference across the nozzle
%rho_g [kg/m3] : gas density
%D_or [m]      : orifice diameter
%
%OUTPUT
%SMD [m]       : Sauter Mean Diameter
%
%EXAMPLE
%[SMD] = nozzle_cone_full2(905,0.0760,0.0323,5*10^5,2.1557,10)

    Delta_p = Delta_p*D_or/sigma_x;
    g       = 9.822;                           %gravitational acceleration [m2/s]
    Oh      = mu_x/(sqrt(D_or*rho_x*sigma_x)); %Ohnesorge number [adim]
    SMD     = D_or*400*Delta_p^-0.75*(rho_g/rho_x/g)^0.2*Oh^(1/7);
end

function [D_g,f_g]= discret_RR(num_cls,delta,beta,alpha)
%This function gives the discretization of the Rosin Rammler volume distribution.
%
%INPUT
%num_cls [adim] : number of droplet classes
%delta [m]      : Rosin Rammler size parameter
%beta [adim]    : Rosin Rammler distribution parameter
%alpha [adim]   : probability that droplet diameters fall below d_min
%                 (tail at alpha/2 of probability) and above d_max (tail at 1-alpha/2)
%
%OUTPUT
%D_g [m]        : mean droplet diameter for each class
%f_g [adim]     : frequency for each class
%
%EXAMPLE
%[D_g,f_g]=discret_RR(num_cls,delta,beta,alpha)

%Discretization parameters for Rosin Rammler
    num    = NaN(1,num_cls);                        % weighted mean numerators pre-allocation
    f_g    = NaN(1,num_cls);                        % class frequencies pre-allocation
    d_min  = delta*log(1/(1-(  alpha/2)))^(1/beta); % value below which droplet diameters from the given distribution fall eps/2 percent of the time
    d_max  = delta*log(1/(1-(1-alpha/2)))^(1/beta); % value above which droplet diameters from the given distribution fall eps/2 percent of the time
    d_i    = linspace(d_min,d_max,num_cls+1);       % class extreme values

    %First class
    num(1) = intgr_RR(d_i(2),delta,beta,'lower');
    f_g(1) = F_RR(d_i(2),delta,beta);

    %Intermediate classes (if present)
    if num_cls >= 3
        num(2:end-1) = intgr_RR(d_i(3:end-1),delta,beta,'lower')-intgr_RR(d_i(2:end-2),delta,beta,'lower');
        f_g(2:end-1) = F_RR(d_i(3:end-1),delta,beta)-F_RR(d_i(2:end-2),delta,beta);
    end

    %Last class
    num(end) = intgr_RR(d_i(end-1),delta,beta,'upper');
    f_g(end) = 1-F_RR(d_i(end-1),delta,beta);
    D_g      = num./f_g; %D_g is calculated according to the theorem of the weighted mean
end

function [cum]= F_RR(d,delta,beta)
%This function calculates the cumulative probability F(d) of the Rosin Rammler
%volume distribution. (Fd) is the volume fraction of drops whose diameters
%are smaller than d.
%
%INPUT
% d [m]       : droplet diameter
% delta [m]   : Rosin Rammler size parameter
% beta [adim] : Rosin Rammler distribution parameter
%
%OUTPUT
%cum [adim]  : cumulative probability
%
%EXAMPLE
%[cum]=F_RR(90*10^-6,200*10^-6,3.5)

    cum = 1-exp(-(d/delta).^beta);
end

function [num]= intgr_RR(d,delta,beta,tail)
%This function evaluates the integral of the function f(d)*d where d is
%the droplet diameter and f(d) is the Rosin Rammler probability density
%function. This calculation is performed using the incomplete gamma
%function.
%
%INPUT
%d [m]         : droplet diameter
%delta [m]     : Rosin Rammler size parameter
%beta [adim]   : Rosin Rammler distribution parameter
%tail [string] : input parameter of the incomplete gamma function
%                'lower' for the lower incomplete gamma function
%                'upper' for the upper incomplete gamma function
%
%OUTPUT
%num [adim]    : value of the integral
%
%EXAMPLE
%[num]=intgr_RR(90*10^-6,150*10-6,3.5,'lower')

    lim = (d/delta).^beta; %upper or lower limit of integration
    a   = 1/beta+1;
    num = delta*gammainc(lim,a,tail)*gamma(a);
end