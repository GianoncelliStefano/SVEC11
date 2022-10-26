function [Q_g_RS_dis,Q_g_PE_dis,Q_g_VS_dis,Q_g_RS_suction,Q_g_PE_suction,Q_g_VS_suction,Q_g_PE_closed_cell,Q_g_VS_closed_cell,Fleak_g,Fleak_g_PE_suc,Fleak_g_VS_suc,Fleak_g_PE_cell_closed,Fleak_g_VS_cell_closed,m_g,q_gas_Id,q_gas,T_l,T_g,T_mix,T_lmean,p_comp,rho_o,nu_o,sigma_o,m_inj_nzl,V_inj_nzl,V_inj_cls,epsm,DU,n_pltr,Dg,fg,cls_indx,M,m_gas_Id,m_gas,m_liq_Id,m_liq,Q_th,R_g,fOK,theta_vane,theta_SucOpen,theta_SucClose,theta_DisOpen,theta_DisClose,Flow_RS,Flow_VS,Flow_PE] = ...
    S3_Thermodynamics(process,D,d,d_hub,L,BI,s,RSclr,VSclr,PEclr,TgAngle,NOZZLES,rho_ref,T_ref,nu_40,nu_100,c_l,c_v,k_g,k_l,MM_g,T_suc,T_0,p_0,p_suc,p_del,mu_g,V_cell,V_comp,pos_SucOpen,pos_SucClose,pos_DisOpen,pos_DisClose,pos_End,Npt_cell,theta,theta_vane,omega,toll_d,Gamma,c,n_van,rpm,model_g,model_full,name4prop,molcomp_g,propSource,MUmodel_RS,MUmodel_VS,MUmodel_PE,RSmodel,VSmodel,PEmodel,IterMax,UndRlx,fLKG,fLKG_in)
% This function calculates temperature and pressure of gas and liquid
% during the compression process, heat transfer,
% specific heat for the polytropic process and the polytropic index.
%
% INPUT
% process             : process performed (1 - compression    2: expansion)
% D [m]               : stator diameter
% d [m]               : rotor diameter
% d_hub [m]           : shaft diameter
% L [m]               : rotor length
% BI [m]            : real vane excursion measured on vane axis
% s [m]             : vane thickness
% RSclr [m]           : rotor-stator clearance width
% VSclr [m]         : vane-end-plate clearance size
% TgAngle [rad]       : rotor-stator tangency angle
% NOZZLES             : structure of nozzles data
% rho_ref [kg/m3]     : reference density for lubricating oil
% T_ref [K]           : reference temperature for rho_ref
% nu_40 [cSt]         : oil kinematic viscosity @40°C
% nu_100 [cSt]        : oil kinematic viscosity @100°C
% c_l [J/kgK]         : specific heat for lubricating oil
% c_v [J/kgK]         : gas specific heat at constant volume
% k_g [W/mK]          : thermal conductivity of gas
% k_l [W/mK]          : thermal conductivity of oil
% MM_g [kg/mol]       : molar mass of gas
% T_suc [K]           : initial temperature of gas
% T_0 [K]             : reference temperature for energy balance
% p_0 [Pa]            : reference pressure
% p_suc [Pa]          : suction pressure
% p_del  [Pa]         : delivery pressure
% mu_g [Pa s]         : gas dynamic viscosity
% V_cell [m3]         : cell volume
% V_comp [m3]         : volumes vector of chamber during compression
% pos_SucOpen [-]     : suction opening discretized angular position in vector theta
% pos_SucClose [-]    : suction closing discretized angular position in vector theta
% pos_DisOpen [-]     : delivery opening discretized angular position in vector theta
% pos_DisClose [-]    : delivery opening discretized angular position in vector theta
% pos_End             : tangency ending discretized angular position in vector theta
% Npt_cell [-]        : number of discretization points of the cell
% theta [rad]         : discretized angulòar position between (-gamma : 2pi)
% theta_vane [rad]    : vane angular positions between (0:2pi)
% omega [rad/s]       : angular velocity
% toll_d [-]          : discretization tolerance
% Gamma [rad]         : cell angular aperture
% c [-]               : geometry index
%   1                     : circular stator
%   2                     : elliptical stator
% n_van [-]           : number of vanes
% rpm [rpm]           : shaft angular speed
% model_g             : thermodynamic model - 'simple' : simplified thermodynamics | 'full': full thermodynamics
% model_full          : thermodynamic model - 'fullTP' : mono-component, single-phase / multi-component, single-phase or multi-phase | 'fullTV' : mono-component, two-phase system
% name4prop           : name of fluid used for real gas thermodynamic property computetion
% molcomp_g           : gas molar composition (expressed as molar fraction of each component)
% propSource          : source for gas properties - 'direct' : ThermoPhysProps | 'REFPROP': link to REFPROP
% MUmodel_RS [string]    : viscosity model for two-phase gas-liquid mixtures
% MUmodel_VS [string] : viscosity model for two-phase gas-liquid mixture
% RSmodel [string]    : leakage model for rotor-stator clearance
% VSmodel [string]    : vane-end-plate leakage model
% IterMax [-]         : Max number of iteration for leakage computation
% UndRlx [-]          : under-relaxation factor for ierative loops
% fLKG [-]            : leakage flag
% fLKG_in [-]         : internal leakages flag
%
% OUTPUT
% m_g [kg]            : mass of gas in the first closed chamber
% q_gas_Id [m3/s]     : ideal gas volumetric flow rate
% q_gas [m3/s]        : real gas volumetric flow rate
% T_l [K]             : liquid temperature during closed chamber phase
% T_g [K]             : gas temperature during closed chamber phase
% T_mix [K]           : adiabatic mixing temperature during closed chamber phase
% T_lmean [K]         : oil adiabatic mixing temperature for each discretization step
% p_comp [Pa]         : gas pressure during closed chamber phase
% rho_o [kg/m3]       : density of injected oil for each nozzles
% nu_o [cSt]          : cinematic viscosity of injected oil for each nozzles
% sigma_o [N/m]       : surface tension of injected oil for each nozzles
% m_inj_nzl [kg]      : mass of oil injected by each nozzle
% V_inj_nzl [m3]      : mass of oil injected by each nozzle
% V_inj_cls [m3]        : volume of oil injected by each nozzle class
% epsm [-]            : ratio between mass of injected oil and gas mass trapped in a cell for each injector
% DU [J]              : variation of internal energy for gas-liquid mixture
% n_pltr [-]          : polytropic index
% Dg [m]              : droplet diameters for each class
% fg [-]              : droplet frequency for each class
% cls_indx [-]        : index of disclosed nozzles classes
% M [-]               : number of injectors
% m_gas_Id [kg/s]     : ideal gas mass flow rate (no leakages)
% m_gas [kg/s]        : mass flow rate of gas
% m_liq_Id [kg/s]     : total liquid mass flow rate
% m_liq [kg/s]        : ideal liquid mass flow rate
% Q_th [W]            : Heat power between gas and liquid
% R_g [J/kgK]         : specific gas constant
% fOK [bool]          : OK flag (1 - no problem occourred, 0 - a problem have been spotted, SVEC will exit the simulation)
%
% NOTE : Hypothesis of adiabatic walls -> DU == Lvano
% HISTORY : V10.1_DeFranco_Genoni_Gianoncelli: review of suction and delivery times (further information available in the thesis) 
    %% DEFINITIONS %%
    fOK         = 1;                                 % process control flag (1- no problem occoured    0- an error eccoured)
    lkg_cls = 1;
    VleakIguess = V_comp(1)*c*rpm*n_van/60*10^(-5);  % firt guess of rotor-stator oil leakage volumetric flow rate [m3/s]
    
    % Extractions from theta (!!! in the future insert also theta_suc e theta_dis for a full 360° revolution)
    theta_comp  = theta(pos_SucClose:pos_DisOpen-Npt_cell)'; % vane angular positions during the compression process [rad]
    theta_i     = theta(pos_SucClose);                       % starting angular position of compression process [rad]
    theta_f     = theta(pos_DisOpen-Npt_cell);               % end angular position of compression process [rad]
    theta_SucOpen = theta_vane(pos_SucOpen-Npt_cell);       % suction open [rad]   
    theta_SucClose   = theta_vane(pos_SucClose-Npt_cell);   % suction close [rad]
    theta_DisOpen = theta_vane(pos_DisOpen-Npt_cell);       % delivery open [rad]
    theta_DisClose    = theta_vane(pos_DisClose-Npt_cell);  % delivery close [rad]
    angle_dis = theta_vane(end)-theta_vane(pos_DisOpen-2*Npt_cell);  % Angle in discharge swept by reference blade [rad] 
    angle_suc = theta_vane(pos_SucClose-Npt_cell);                   % Angle in suction swept by reference blade [rad]
   
    % Process parameters
    Npt_comp        = length(V_comp);                        % discretization points of closed cell process [-]
    t_comp          = (theta_f-theta_i)/omega;               % duration of the compression process [s]
    dt              = t_comp/Npt_comp;                       % time step [s]
    dt_vano         = Gamma/omega;                           % time needed for a vane to sweep a cell [s]
    dt_angle_dis  = angle_dis/omega;                         % A closed cell see discharge port for this amount of time [s]
    dt_angle_suc  = angle_suc/omega + dt_vano;               % A closed cell see suction port for this amount of time(include both phase I and phase II in suction, see S3_LeakagePE) [s]
  
    % Thermodynamics
    R_g        = SX_Constant({'UniGasConstant'})/MM_g;  % specific gas constant [J/ kg K]
    c_p        = c_v + R_g;                             % gas specific heat at costant pressure [J/kgK]
    v_spec_0   = R_g*T_0/p_0;                           % specific volume of reference [m^3/kg]
    clear theta_f t_comp

    %% LEAKAGES PREAMBLE %%
    % COMPRESSORS: leakage enabled & oil injection.
    % In Compressors, leaked oil form rotor-stator clearance is obtained
    % via a fictious nozzle, that account for the termal effect of the leaking oil
    if fLKG && process == 1 && sum(NOZZLES.f_nz)>=1  
        VleakRS_E              = 0;     % oil flow rate from suction to delivery via Rotor-stator clearance [m3/s](only expander, but needed for overall calculation)
        NOZZLES.f_nz           = vertcat(1,NOZZLES.f_nz);
        NOZZLES.theta_nz       = vertcat(0,NOZZLES.theta_nz);
        NOZZLES.num_nz         = vertcat(1,NOZZLES.num_nz);
        NOZZLES.p_inj          = vertcat(NaN,NOZZLES.p_inj);
        NOZZLES.Tl_in          = vertcat(mean(NOZZLES.Tl_in),NOZZLES.Tl_in);
        NOZZLES.num_cls        = vertcat(1,NOZZLES.num_cls);
        NOZZLES.alpha_nz       = vertcat(NaN,NOZZLES.alpha_nz);
        NOZZLES.type_nz        = vertcat('kf',NOZZLES.type_nz);
        NOZZLES.name_nz        = vertcat('ud',NOZZLES.name_nz);
        NOZZLES.specs          = horzcat({struct},NOZZLES.specs);
        NOZZLES.specs{1}.D_g   = 0.001;
        NOZZLES.specs{1}.f_g   = 1;
        NOZZLES.specs{1}.V_inj = VleakIguess; % first guess of compressos 'injected' rotor-stator leakage volumetric flow rate [m3/s]
    % EXPANDER: Leakage enabled & oil injection.
    % in expander, leakage is obtained by recalculating the initial oil
    % volume after leaking oil has been calculated
    elseif fLKG && process == 2 && sum(NOZZLES.f_nz)>=1 
        VleakRS_E = VleakIguess;  % first guess of expander rotor-stator leakage volumetric flow rate [m3/s]
    % No oil injection or no leakages
    else
        VleakRS_E = 0;
        lkg_cls   = [];
    end

    % Unpack basic nozzle data
    Nu       = NOZZLES.Nu;       % Nusselt number [-]
    num_cls  = NOZZLES.num_cls;  % number of classes per each nozzle
    Tl_in    = NOZZLES.Tl_in;    % nozzles injection temperature [K]
    theta_nz = NOZZLES.theta_nz; % nozzles angular position [rad] 
    
    %% PRELIMINARY CALCULCATION %%
    M          = length(theta_nz);                     % number of nozzles
    K          = sum(num_cls);                         % total number of classes
    [~,seq_nz] = sort(theta_nz);                       % indexes of nozzles as if they were sorted
    cls_indx   = repelem(seq_nz,num_cls);              % index of disclosed nozzles classes

    % Nozzle positions in the compression process
    cnd_suc    = (theta_i+Gamma/2)>theta_nz;                                                % condition to have a noozle during suction
    nzls_suc   = sum(cnd_suc);                                                              % number of active nozzles present before suction close plus fictious rotor-stator clearance nozzle
    [~,pos_nz] = arrayfun(@(x) min(abs(theta_comp-(theta_nz(x)-Gamma/2))), (nzls_suc+1):M); % indexes of theta_comp array where a new nozzles is disclosed, during closed cell phase
    clear seq_nz theta_nz_deg theta_comp

    %% SUCTION %%
    Loop_Out  = 1; % outer iteration flag
    Niter_Out = 1; %Iteration counter

    while Loop_Out
        % preallocation
        V_l        = NaN(Npt_comp,1);  % volume of injected oil for each discretization point [m3]
        m_l        = NaN(Npt_comp,1);  % mass of injected oil for each discretization point [kg]
        V_inj_nzl  = NaN(1,M);         % volume of oil injected by each nozzle  [m3]
        m_inj_nzl  = NaN(1,M);         % mass of oil injected by each nozzle [kg]
        V_inj_cls  = NaN(1,K);         % volume of injected oil for each class [m3]
        m_inj_cls  = NaN(1,K);         % mass of injected oil for each class [kg]
        V_cell_cls = NaN(Npt_comp,K);  % volume of in the cell for each class for each class discretization step [m3]
        m_cell_cls = NaN(Npt_comp,K);  % mass of in the cell oil for each class for each class discretization step [kg]
        V_g        = NaN(Npt_comp,1);  % volume of gas for each discretization point [m3]
        m_g        = NaN(Npt_comp,1);  % mass of gas for each discretization point [kg]
        v_spec     = ones(Npt_comp,1); % sepcific volume for each discretization point [m^3/kg]
        quality    = NaN(Npt_comp,1);  % gas quality for each discretization point [kg_v/kg_tot]
        T_l        = NaN(Npt_comp,K);  % liquid temperature for each discretization point 
        T_g        = NaN(Npt_comp,1);  % gas temperature for each discretization poiand each class [K]
                                       %  rows=iteration step - columns=class
        p_comp     = NaN(Npt_comp,1);  % gas pressure for each discretization point [Pa]
        Pind_ist   = NaN(Npt_comp,1);  % indicated power for each discretization point [W]
        Z          = NaN(Npt_comp,1);  % compressibility factor [-]
        cx         = NaN(Npt_comp,1);  % specific heat of politropic process [J/kgK]
        n_pltr     = NaN(Npt_comp,1);  % politropic index
        A          = NaN(K+1,K+1);     % global coefficient matrix A
        B          = NaN(K+1,1);       % global constant terms B
        dUgas      = NaN(Npt_comp-1,1);  % Infinitesimal variation of gas internal energy [J]
        dUliq      = NaN(Npt_comp-1,1);  % Infinitesimal variation of oil internal energy [J]
        dH_g       = NaN(Npt_comp-1,1);  % Infinitesimal variation of gas enthalpy energy due to leakages [J]
        dH_o       = NaN(Npt_comp-1,1);  % Infinitesimal variation of oil enthalpy energy due to leakages [J]
        L_g        = NaN(Npt_comp-1,1);
        dQ         = NaN(Npt_comp,K);
        C_l        = NaN(Npt_comp,K);
        Dg         = NaN(1,K);         % droplet diameters for each class [m]
        fg         = NaN(1,K);         % class frequencies [adim]
        map_cls    = zeros(1,K);       % class mapping vector
        rho_o      = NaN(1,M);         % oil density for each nozzle [kg/m3]
        nu_o       = NaN(1,M);         % oil kinematic viscosity for each nozzle [cSt]
        mu_o       = NaN(1,M);         % oil dynamic viscosity for each nozzle [Pa s]
        sigma_o    = NaN(1,M);         % oil surface tension for each nozzle [N m]

        % conductivity and compressibility factor during suction
        switch model_g
            case 'simple'
                k_suc = k_g;
                Z_suc = 1;
            case 'full'
                [quality(1)] = ThermoPhysProps('Q',T_suc,p_suc,name4prop,molcomp_g,propSource,0,'PR','ChungEtal',1e-3,80);
                if quality(1) == 0
                    warning('S3_Thermodynamics:quality','Two phase inlet flow: undefined system with T_suc and P_suc (correlated).');
                    SX_Logfile ('e',{lastwarn});
                    fOK = 0;
                else
                    [k_suc,Z_suc] = ThermoPhysProps('kv,zf',T_suc,p_suc,name4prop,molcomp_g,propSource,0,'PR','ChungEtal',1e-3,80);
                end
        end
        
        % suction preallocation
        for j = 1:(nzls_suc)
            if fOK
                map_cls(cls_indx==j) = 1;        % increment counter for total disclosed nozzle classes
                T_l(1,cls_indx==j)   = Tl_in(j); % update injection temperature matrix
                [rho_o(j),nu_o(j),mu_o(j),sigma_o(j),fOK] = S3_OilProperties(Tl_in(j),T_ref,rho_ref,nu_40,nu_100); % unpack oil properties for current nozzle
                
                if strcmp(NOZZLES.type_nz(j), 'un')
                    m_g_ideal    = p_suc*V_comp(1)/(R_g*T_suc*Z_suc);               % ideal trapped gas mass [kg]
                    D_g          = NOZZLES.specs{j}.D_g;                            % class droplet diameter, from user input [m]
                    f_g          = NOZZLES.specs{j}.f_g;                            % class frequency [-]
                    m_inj_nzl(j) = NOZZLES.specs{j}.eps_m*m_g_ideal; % mass of oil injected in cell [kg]
                    V_inj_nzl(j) = m_inj_nzl(j)/rho_o(j);                           % Volume of injected oil for curent nozzle [m3]
                else
                    [V_nz,m_nz,~,D_g,f_g,fOK] = S3_Nozzles(NOZZLES,j,rho_o(j),nu_o(j),mu_o(j),sigma_o(j),p_suc,T_suc,R_g);
                    V_inj_nzl(j)              = V_nz*dt_vano;  % Volume of injected oil for curent nozzle [m3]
                    m_inj_nzl(j)              = m_nz*dt_vano;  % Mass of injected oil for curent nozzle [kg]
                end
            end
           
            % Update global nozzle index properties
            Dg(cls_indx==j)        = D_g;              % update array Dg with value from disclosed nozzle
            fg(cls_indx==j)        = f_g;              % update array fg with value from disclosed nozzle
            V_inj_cls(cls_indx==j) = V_inj_nzl(j)*f_g; % update total injected oil volume per class with value of disclosed nozzle
            m_inj_cls(cls_indx==j) = m_inj_nzl(j)*f_g; % update global injected oil mass per class with value of disclosed nozzle
        end
        
        map_cls_suc = map_cls; % counter for total disclosed nozzle classes in suction
        
        % ===============================================================
        % ===============================================================
        % PREALLOCATION
        % Thermodynamics
        p_comp_C        = p_suc*ones(Npt_comp,1);                   % pressure vector at previous inner iteration [Pa]  p_comp_C = p_suc.*(V_comp(1)./V_comp').^(c_p/c_v); (option: considering an adiabatic process as first guess, but problems occur in case of oil free)
        T_leak_g        = T_suc*ones(Npt_comp+2*Npt_cell,1);        % Temperature vector for gas leaking enthalpy flow along the compression [K]
        T_leak_o        = zeros(Npt_comp+2*Npt_cell,1);             % Temperature vector for oil leaking enthalpy flow along the compression [K]
        p_leak          = p_suc*ones(Npt_comp+2*Npt_cell,1);        % Pressure vector for leaking enthalpy flow along the compression [Pa]
        v_spec_leak     = ones(Npt_comp+2*Npt_cell,1);               % Specific volume vector for leaking enthalpy flow along the compression [m^3/kg]
        rho_o_Lkg23     = ones(1,1);                                % Leaking oil density [kg/m^3]
        % Mass flow rate [kg/s]
        Fleak_g=zeros(1,1); 
        Fleak_g_PE_suc      = zeros(1,1);                   % leaking gas mass flow rate in suction through path 2 [kg/s]
        Fleak_g_VS_suc      = zeros(1,1);                   % leaking gas mass flow rate in suction through path 3 [kg/s]
        Fleak_g_VS_in_TV    = zeros(Npt_comp-1,1);          % leaking gas mass flow rate TO the cell through path 3 through trailing vane[kg/s]
        Fleak_g_VS_in_LV    = zeros(Npt_comp-1,1);          % leaking gas mass flow rate TO the cell through path 3 through leading vane[kg/s]
        Fleak_g_VS_out      = zeros(Npt_comp-1,1);          % leaking gas mass flow rate FROM the cell through path 3 [kg/s]
        Fleak_g_VS_dis      = zeros(1,1);                   % leaking gas mass flow rate in discharge through path 3 [kg/s]
        Fleak_g_PE_dis      = zeros(1,1);                   % leaking gas mass flow rate in discharge through path 2 [kg/s]
        Fleak_g_PE_net      = zeros(Npt_comp-1,1);          % net leaking gas mass flow rate in discharge through path 2 [kg/s]
        Fleak_o_VS_in_TV    = zeros(Npt_comp-1,1);          % leaking oil mass flow rate TO the cell through path 3 through trailing vane[kg/s]
        Fleak_o_VS_in_LV    = zeros(Npt_comp-1,1);          % leaking oil mass flow rate TO the cell through path 3 through leading vane[kg/s]
        Fleak_o_VS_out      = zeros(Npt_comp-1,1);          % leaking oil mass flow rate FROM the cell through path 3 [kg/s]
        Fleak_o_VS_suc      = zeros(1,1);
        Fleak_o_VS_dis      = zeros(1,1);                   % leaking oil mass flow rate in discharge through path 3 [kg/s]
        Fleak_o_PE_dis      = zeros(1,1);                   % leaking oil mass flow rate in discharge through path 2 [kg/s]
        Fleak_o_PE_net      = zeros(Npt_comp-1,1);          % net leaking oil mass flow rate along compression process through path 2 [kg/s]
        Fleak_o_net         = zeros(Npt_comp-1,1);          % net leaking oil mass flow rate  through path 2 and 3 [kg/s]
        Fleak_o_cls0_C      = zeros(Npt_comp-1,1);
        Fleak_o_cls0_E      = zeros(Npt_comp-1,1);
        Fleak_g_PE_cell_closed = zeros(1,1);
        Fleak_g_VS_cell_closed = zeros(1,1);     
        Flow_PE = zeros(length(theta_vane),1); 
        Flow_VS = zeros(length(theta_vane),1); 
        Flow_RS = zeros(length(theta_vane),1); 
        Q_g_PE_suction = zeros(1,1);
        Q_g_VS_suction = zeros(1,1);
        Q_g_PE_closed_cell = zeros(1,1);
        Q_g_VS_closed_cell = zeros(1,1);
        
        % Volume flow rate [m^3/s]
        Vleak_o_suc_C       = zeros(1,1);                   % leaking oil volume in suction in compressor [m^3/s]
        Vleak_o_suc_E       = zeros(1,1);                   % leaking oil volume in suction in expander [m^3/s]
        Vleak_o_net         = zeros(Npt_comp-1,1);          % net leaking oil volume in compression [m^3/s]
        Vleak_o_VS_out      = zeros(Npt_comp-1,1);          % leaking oil FROM the cell [m^3/s]
        Vleak_o_cls0_C      = zeros(Npt_comp-1,1);          % leaking oil TO the cell accounting only for oil class 0 [m^3/s]
        Vleak_o_cls0_E      = zeros(Npt_comp-1,1);          % leaking oil TO the cell accounting only for oil class 0 [m^3/s]
        % Mass [kg]
        Dmg_VS_comp         = zeros(Npt_comp,1);          % cumulated gas mass variation due to leakage through path 3 along compression process [kg]
        Dmg_PE_comp         = zeros(Npt_comp,1);          % cumulated gas mass variation due to leakage through path 2 along compression process [kg]
        
        Loop_In  = 1; % outer iteration flag
        Niter_In = 1; %Iteration counter
        % ===============================================================
        % ===============================================================
        while Loop_In
            % index inizialization
            map_cls = map_cls_suc;
            
            % update initial cell oil volume and mass in case of leakage in expander
            V_cell_cls(1,map_cls==1) = V_inj_cls(map_cls==1)*(1-(Vleak_o_suc_E + VleakRS_E)*dt_vano/sum(V_inj_nzl(cnd_suc)));      % initial oil volume in cell for each class [m3]
            % Accounting for leaking oil in suction (all oil leakage is considered to belong to the class introduced with the fictitious nozzle D_g = 1mm, see Franzetti-Persico thesis)
            V_cell_cls(1,lkg_cls) = V_cell_cls(1,lkg_cls) + Vleak_o_suc_C*dt_vano;
            m_cell_cls(1,map_cls==1) = V_cell_cls(1,map_cls==1).*repelem(rho_o(1:nzls_suc),num_cls(1:nzls_suc)); % initial oil mass in cell for each class [kg]
            
            h     = Nu*k_suc./Dg;                          % update heat transfer coefficient for each class [W/m^2K]
            Bi    = h/k_l.*Dg/6;                           % update biot number for each class [adim]
            Al    = pi.*Dg.^2;                             % update droplet surface for each class [m2]
            Vdrop = pi.*Dg.^3/6;                           % update droplet volume for each class [m3]
            % NS    = V_cell_cls(1,1:end)./Vdrop;  update number of droplets for each class [adim]
            
            % air mass calculation
            V_l(1) = sum(V_cell_cls(1,map_cls==1));  % volume of oil trapped in first cell of the machine [m3]
            m_l(1) = sum(m_cell_cls(1,map_cls==1));  % mass of oil trapped in first cell of the machine [kg]
            V_g(1) = (V_comp(1)-V_l(1));             % volume of gas trapped in the first cell of the machine [m3]
            m_g(1) = p_suc*V_g(1)/(R_g*T_suc*Z_suc); % mass of air elaborated by first cell of the machine [kg]
            

            m_g = m_g(1) + Dmg_VS_comp + Dmg_PE_comp;
          
            clear V_nz_tot m_nz_tot D_g f_g alpha num_nz
            
            %% CLOSED-CHAMBER ITERATIVE SOLUTION %%
            % Initialization
            T_g(1)    = T_suc;                                % starting temperature [K]
            p_comp(1) = p_suc;                                % starting pressure [Pa]
            Z(1)      = Z_suc;                                % compressibility factor during suction [-]
            v_spec(1) = V_g(1)/m_g(1);                        % starting specific volume [m^3/kg]
            
            % Looping over discretization points
            for i = 1:(Npt_comp-1)
                if fOK
                    DV = Vleak_o_net(i)*dt;  % inizialization of oil volume injected at current step
                    Dm = Fleak_o_net(i)*dt;  % inizialization of oil mass injected at current step
                    
                    % nozzle evaluation
                    [val,pos] = min(abs(i-pos_nz));
                    if val == 0                          % new nozzle disclosed
                        j = nzls_suc + pos;              % number of the disclosed nozzle
                        map_cls(cls_indx==j) = 1;        % array counter for total disclosed nozzle classes
                        T_l(i,cls_indx==j)   = Tl_in(j); % update injection temperature matrix
                        [rho_o(j),nu_o(j),mu_o(j),sigma_o(j),fOK] = S3_OilProperties(Tl_in(j),T_ref,rho_ref,nu_40,nu_100); % unpack oil properties for current nozzle
                        
                        % extract nozzle-specific data for that nozzle
                        if strcmp(NOZZLES.type_nz(j),'un')
                            m_inj_nzl(j) = NOZZLES.specs{j}.eps_m*m_g(i+1);
                            V_inj_nzl(j) = m_inj_nzl(j)/rho_o(j);
                            D_g          = NOZZLES.specs{j}.D_g;
                            f_g          = NOZZLES.specs{j}.f_g;
                        else
                            [V_nz,m_nz,~,D_g,f_g,fOK] = S3_Nozzles(NOZZLES,j,rho_o(j),nu_o(j),mu_o(j),sigma_o(j),p_comp(i),T_g(i),R_g);
                            V_inj_nzl(j)              = V_nz*dt_vano;  % Volume of injected oil for curent nozzle [m3]
                            m_inj_nzl(j)              = m_nz*dt_vano;  % Mass of injected oil for curent nozzle [kg]
                        end
                        
                        % Update global nozzle index properties
                        DV                     = DV + V_inj_nzl(j);     % volume of oil injected at current step [m3]
                        Dm                     = Dm + m_inj_nzl(j);     % mass of oil injected at current step [kg]
                        Dg(cls_indx==j)        = D_g;                   % update array Dg with value from disclosed nozzle [m]
                        fg(cls_indx==j)        = f_g;                   % update array fg with value from disclosed nozzle [-]
                        V_inj_cls(cls_indx==j) = V_inj_nzl(j)*f_g;      % update total injected oil volume per class with value of disclosed nozzle [m3]
                        m_inj_cls(cls_indx==j) = m_inj_nzl(j)*f_g;      % update global injected oil mass per class with value of disclosed nozzle [kg]
                        
                        % add to cell oil class volume and mass of just disclosed nozzle
                        V_cell_cls(i,cls_indx==j) = V_inj_cls(cls_indx==j);  % oil volume in cell for each class in current time step [m3]
                        m_cell_cls(i,cls_indx==j) = m_inj_cls(cls_indx==j);  % oil mass in cell for each class in current time step [kg]
                        
                        Al    = pi.*Dg.^2;         % update droplet surface for each class [m2]
                        Vdrop = pi.*Dg.^3/6;       % update droplet volume for each class [m3]
                        
                        clear V_nz_tot m_nz_tot D_g f_g
                    end
                    
                    % Volume and mass for each class of nozzle are updated.
                    V_cell_cls(i+1,map_cls==1) = V_cell_cls(i,map_cls==1)*(1+(Vleak_o_cls0_E(i)-Vleak_o_VS_out(i))*dt/sum(V_cell_cls(i,map_cls==1)));                      %-VleakPE +\- VleakVS;   oil volume in cell for each class in current time step [m3]
                    V_cell_cls(i+1,lkg_cls)    = V_cell_cls(i+1,lkg_cls) + Vleak_o_cls0_C(i)*dt;                                                          %-VleakPE +\- VleakVS;   oil volume in cell for each class in current time step [m3] accounting for class 0
                    m_cell_cls(i+1,map_cls==1) = m_cell_cls(i,map_cls==1)*(1+(Fleak_o_cls0_E(i)-Fleak_o_VS_out(i))*dt/sum(m_cell_cls(i,map_cls==1)));                      %-MleakPE +\- MleakVS;   oil mass in cell for each class in current time step [kg]
                    m_cell_cls(i+1,lkg_cls)    = m_cell_cls(i+1,lkg_cls) + Fleak_o_cls0_C(i)*dt;                                  %-VleakPE +\- VleakVS;   oil mass in cell for each class in current time step [kg] accounting for class 0
                    
                    % update values
                    V_l(i+1) = V_l(i)+ DV;
                    m_l(i+1) = m_l(i)+ Dm;
                    V_g(i+1) = V_comp(i+1)-V_l(i+1);
                    v_spec(i+1) = V_g(i+1)/m_g(i+1);
                    
                    % =====================================================
                    % =====================================================
                    
                    NS    = V_cell_cls(i+1,1:end)./Vdrop;  % update number of droplets for each class [adim]
                    
                    % conductivity factor at current step
                    switch model_g
                        case 'simple'
                            k_real = k_g;
                        case 'full'
                            switch model_full
                                case 'fullTP'
                                    [k_real] = ThermoPhysProps('kv',T_g(i),p_comp(i),name4prop,molcomp_g,propSource,0,'PR','ChungEtal',1e-3,80);
                                case 'fullTV'
                                    [quality(i),rho_liq,rho_vap] = ThermoPhysProps('Q,rhol,rhov',T_g(i),v_spec(i),name4prop,molcomp_g,'REFPROP_TV',0);
                                    % specific volumes of saturated liquid and vapour phase: corrective coefficient are defined inorder to avoid numerical problem with
                                    % Thermophysprops solution for saturated phase:
                                    fcorr = 0.01;
                                    v_liq = 1/rho_liq*(1-fcorr);
                                    v_vap = 1/rho_vap*(1+fcorr);
                                    
                                    if quality(i) > 0 && quality(i) < 1
                                        [k_liq] = ThermoPhysProps('kl',T_g(i),v_liq,name4prop,molcomp_g,'REFPROP_TV',0);
                                        [k_vap] = ThermoPhysProps('kv',T_g(i),v_vap,name4prop,molcomp_g,'REFPROP_TV',0);
                                        % two-phase conductivity (EMT method by Awad et al. - 2008, see Franzetti-Persico thesys
                                        k_real = 1/4*(sum([(3*quality(i)-1).*k_vap,(3*(1-quality(i))-1).*k_liq],2,'omitnan') + sqrt(sum([sum([(3*quality(i)-1).*k_vap,(3*(1-quality(i))-1).*k_liq],2,'omitnan').^2,8*k_liq.*k_vap],2,'omitnan')));
                                    elseif quality(i) >= 1
                                        [k_real] = ThermoPhysProps('kv',T_g(i),v_spec(i),name4prop,molcomp_g,'REFPROP_TV',0);
                                    else
                                        [k_real] = ThermoPhysProps('kl',T_g(i),v_spec(i),name4prop,molcomp_g,'REFPROP_TV',0);
                                        warning('S3_Thermodynamics:conductivity','The flow is fully condensed.');
                                        SX_Logfile ('w',{strrep(lastwarn,'%','%%')});
                                    end
                            end
                    end
                    
                    h     = Nu*k_real./Dg;     % update heat transfer coefficient for each class [W/m^2K]
                    Bi    = h/k_l.*Dg/6;       % update biot number for each class [adim]
                    
                    switch model_g
                        case 'simple'
                            % Global coefficient matrix A is updated
                            A(1,1)         = m_g(i+1)*c_v;
                            A(1,2:end)     = m_cell_cls(i+1,1:end).*c_l;
                            A(2:end,1)     = -h.*Al.*NS.*dt;
                            A(2:end,2:end) = diag(m_cell_cls(i+1,1:end).*c_l+h.*Al.*NS.*dt);
                            
                            % Step matrix A1 and B1 are generated
                            A1 = A(1:sum(map_cls)+1,1:sum(map_cls)+1); % A1 is generated from A
                            B1 = B(1:sum(map_cls)+1,1);                % B1 is initialized to a NaN vector of proper dimension
                            
                            % Extract pertinet coefficient matrix and constant terms
                            k         = m_cell_cls(i+1,map_cls==1).*c_l.*T_l(i,map_cls==1);
                            bg_1      = m_g(i+1)*c_v*T_g(i) - p_comp(i)*(V_g(i+1)-V_g(i)) + sum(k(1:sum(map_cls)));
                            bg_2      = (Fleak_g_VS_in_LV(i)*(T_leak_g(i+1+2*Npt_cell)-T_0) + Fleak_g_VS_in_TV(i)*(T_leak_g(i+1)-T_0) - Fleak_g_VS_out(i)*(T_leak_g(i+1+Npt_cell)-T_0))*c_p*dt + Fleak_g_PE_net(i)*c_p*dt*(mean(T_leak_g(1:2000:end))-T_0);                    % (mean(T_leak_g)-T_0);
                            bg_3      = (Fleak_o_VS_in_LV(i)*(T_leak_o(i+1+2*Npt_cell)-T_0) + Fleak_o_VS_in_TV(i)*(T_leak_o(i+1)-T_0) - Fleak_o_VS_out(i)*(T_leak_o(i+1+Npt_cell)-T_0))*c_l*dt + Fleak_o_PE_net(i)*c_l*dt*(mean(T_leak_o(1:2000:end))-T_0);                  % (mean(T_leak_o)-T_0);
                            bg_4      = p_leak(i+1+2*Npt_cell)*Fleak_o_VS_in_LV(i)*dt/rho_o_Lkg23 + p_leak(i+1)*Fleak_o_VS_in_TV(i)*dt/rho_o_Lkg23 - p_leak(i+1+Npt_cell)*Fleak_o_VS_out(i)*dt/rho_o_Lkg23 + mean(p_leak(1:2000:end))*Fleak_o_PE_net(i)*dt/rho_o_Lkg23;           % mean(p_leak)
                            
                            % B1 is created
                            B1(1)     = bg_1 + bg_2 + bg_3 + bg_4;
                            B1(2:end) = k(1:sum(map_cls));
                            
                            % solve linear system by matrix divide
                            T = A1\B1;
                            
                        case 'full'
                            % vectors definition to semplify fsolve call:
                            % (1): leakage (through vane side) TO the cell via leading vane
                            % (2): leakage (through vane side) TO the cell via trailing vane
                            % (3): leakage (through vane side) FROM the cell
                            % (4): net leakage (through rotor-end-plate)
                            Fleak_g_23      = [Fleak_g_VS_in_LV(i);           Fleak_g_VS_in_TV(i);     Fleak_g_VS_out(i);           Fleak_g_PE_net(i)];
                            Fleak_o_23      = [Fleak_o_VS_in_LV(i);           Fleak_o_VS_in_TV(i);     Fleak_o_VS_out(i);           Fleak_o_PE_net(i)];
                            Tleak_g_23      = [T_leak_g(i+1+2*Npt_cell);      T_leak_g(i+1);           T_leak_g(i+1+Npt_cell);      mean(T_leak_g(1:2000:end))];
                            Tleak_o_23      = [T_leak_o(i+1+2*Npt_cell);      T_leak_o(i+1);           T_leak_o(i+1+Npt_cell);      mean(T_leak_o(1:2000:end))];
                            p_leak_23       = [p_leak(i+1+2*Npt_cell);        p_leak(i+1);             p_leak(i+1+Npt_cell);        mean(p_leak(1:2000:end))];
                            v_spec_leak_23  = [v_spec_leak(i+1+2*Npt_cell);   v_spec_leak(i+1);        v_spec_leak(i+1+Npt_cell);   mean(v_spec_leak(1:2000:end))];
                            
                            % solve non-linear system through FSOLVE
                            Ki     = sum(map_cls); % liquid classes at current step
                            NotNaN = ~isnan(T_l(i,:));
                            % initial guess of temperatures (same as step before)
                            T_guess         = NaN(1,Ki+1);
                            T_guess(1)      = T_g(i);
                            T_guess(2:Ki+1) = T_l(i,NotNaN);
                            options.Display = 'off';
%                             options.MaxFunctionEvaluations = 1e5;
%                             options.ConstraintTolerance = 1e-9;
%                             options.FunctionTolerance = 1e-9;
%                             options.StepTolerance     = 1e-9;
                            T  = fsolve(@(T) system(T,Ki,T_g(i),T_l(i,NotNaN),Tleak_g_23',Tleak_o_23',T_0,p_0,v_spec_0,p_comp(i),p_leak_23',rho_o_Lkg23,V_g(i:i+1),...
                                m_g(i+1),v_spec(i:i+1),v_spec_leak_23',m_cell_cls(i+1,NotNaN),Fleak_g_23',Fleak_o_23',c_l,h(NotNaN),Al(NotNaN),...
                                NS(NotNaN),dt,model_full,name4prop,molcomp_g,propSource),T_guess,options);
                            clear Ki T_guess NotNaN
                    end
                    % Extract solution
                    T_g(i+1)              = T(1);
                    T_l(i+1,(map_cls==1)) = T(2:end);
                    % compressibility factor at current step
                    switch model_g
                        case 'simple'
                            Z(i+1)      = 1;
                        case 'full'
                            switch model_full
                                case 'fullTP'
                                    [Z(i+1)] = ThermoPhysProps('zf',T_g(i+1),p_comp(i),name4prop,molcomp_g,propSource,0,'PR','ChungEtal',1e-3,80);
                                case 'fullTV'
                                    [Z(i+1)] = ThermoPhysProps('zf',T_g(i+1),v_spec(i+1),name4prop,molcomp_g,'REFPROP_TV',0);
                            end
                    end
                    
                    p_comp(i+1) = m_g(i+1)*R_g*T_g(i+1)/V_g(i+1).*Z(i+1);
                    
                    % update heat from oil
                    dQ(i,:)  = h.*Al.*NS.*map_cls.*(T_g(i)-T_l(i,:)); % heat power added/subtracted to gas due to injected oil, for each class of disclosed noozles
                    
                    % =====================================================
                    % =====================================================
                    switch model_g
                        case 'simple'
                            dUgas(i) = m_g(i+1)*c_v*(T_g(i+1)-T_g(i));
                            dH_g(i)  = bg_2;        % Gas enthalpy variation along the compression control volume due to leakages through path 2 and 3 [J]
                        case 'full'
                            switch model_full
                                case 'fullTP'
                                    [u_g]    = ThermoPhysProps('u',[T_g(i),T_g(i+1)],[p_comp(i),p_comp(i)],name4prop,molcomp_g,propSource,0,'PR','ChungEtal',1e-3,80);
                                    [h_g]    = ThermoPhysProps('h',[Tleak_g_23',T_0],[p_leak_23',p_0],name4prop,molcomp_g,propSource,0,'PR','ChungEtal',1e-3,80);
                                case 'fullTV'
                                    [u_g]    = ThermoPhysProps('u',[T_g(i),T_g(i+1)],[v_spec(i),v_spec(i+1)],name4prop,molcomp_g,'REFPROP_TV',0);
                                    [h_g]    = ThermoPhysProps('h',[Tleak_g_23',T_0],[v_spec_leak_23',v_spec_0],name4prop,molcomp_g,'REFPROP_TV',0);
                            end
                            
                            dUgas(i) = m_g(i+1)*diff(u_g);
                            dH_g(i)  = (Fleak_g_23(1)*(h_g(1)-h_g(5)) + Fleak_g_23(2)*(h_g(2)-h_g(5)) - Fleak_g_23(3)*(h_g(3)-h_g(5)))*dt + Fleak_g_23(4)*dt*(h_g(4)-h_g(5));   
                    end
                    
                    dUliq(i) = sum(m_cell_cls(i+1,1:end).*c_l.*(T_l(i+1,:)-T_l(i,:)),'omitnan');
                    dH_o(i)  = (Fleak_o_VS_in_LV(i)*(T_leak_o(i+1+2*Npt_cell)-T_0) + Fleak_o_VS_in_TV(i)*(T_leak_o(i+1)-T_0) - Fleak_o_VS_out(i)*(T_leak_o(i+1+Npt_cell)-T_0))*c_l*dt + Fleak_o_PE_net(i)*c_l*dt*(mean(T_leak_o(1:2000:end))-T_0)+...
                               + p_leak(i+1+2*Npt_cell)*Fleak_o_VS_in_LV(i)*dt/rho_o_Lkg23 + p_leak(i+1)*Fleak_o_VS_in_TV(i)*dt/rho_o_Lkg23 - p_leak(i+1+Npt_cell)*Fleak_o_VS_out(i)*dt/rho_o_Lkg23 + mean(p_leak(1:2000:end))*Fleak_o_PE_net(i)*dt/rho_o_Lkg23;
                    L_g(i)   = - p_comp(i)*(V_g(i+1)-V_g(i));   % energy exchanged between machine and cell
                    % =====================================================
                    % =====================================================
%                   
                    C_l(i,:) = m_cell_cls(i,1:end).*map_cls.*c_l;                 % heat Capacity for each class disclosed noozles
                end
            end
            clear k T k_real
            % update heat from oil for last step (i=Npt)
            if strcmp(model_g, 'full')
                if strcmp(model_full, 'fullTV')
                    quality(end) =  ThermoPhysProps('Q',T_g(i),v_spec(i),name4prop,molcomp_g,'REFPROP_TV',0);
                end
            end
            NS    = V_cell_cls(end,1:end)./Vdrop;  % update number of droplets for each class [adim]
            dQ(end,:) = h.*Al.*NS.*(T_g(end)-T_l(end,:));
            C_l(end,:)= m_cell_cls(end,1:end).*c_l;
            
            %% ADIABATIC MIXING TEMPERATURE %%
            Q_l     = (T_l .* C_l);                                    % Oil entalphy for each discretization step for each class [J]
            Q_lmean = sum(Q_l,2,'omitnan');                            % Mean oil entalphy for each discretization step [J]
            C_lmean = sum(C_l,2,'omitnan');                            % Mean oil heat capacity for each discretization step [J/K]
            T_lmean = Q_lmean./C_lmean;                                % Oil adiabatic mixing temperature for each discretization step
            T_mix   = (Q_lmean + T_g.*m_g.*c_v)./(C_lmean + m_g*c_v);  % Gas+oil adiabatic mixing temperature
            
            %% INTERNAL LEAKAGES %%
            % leakages through clearance between rotor-end-plates (path 2) and vane-sides (path 3)
            if fLKG_in
                p_new = p_comp; % pressure along compression, new iteration [Pa]
                
                % Compute error between results of current iteration respect to previous one
                errp_comp = max(abs((p_new-p_comp_C)./p_comp_C)); % pressure along compression - error on inner iteration (VS leakage) [-]
                
                % Check iteration error, if it is low, exit outer cycle otherwise update leaking oil mass and density and re-iterate
                if errp_comp <= toll_d
                    Loop_In = 0;
                else
                    p_comp_C = UndRlx*p_new + (1-UndRlx)*p_comp_C; % update pressure vector
                    
                    % Temperature and pressure vector for gas leaking enthalpy flow along the compression
                    T_leak_g = [T_g(1)*ones(Npt_cell,1); T_g; T_g(end)*ones(Npt_cell,1)];
                    
                    if isnan(p_del)                                           % definition of delivery pressure
                        p_out_C = p_comp_C(end);
                    else
                        p_out_C = p_del;
                    end
                    
                    p_leak      = [p_comp_C(1)*ones(Npt_cell,1); p_comp_C; p_out_C(end)*ones(Npt_cell,1)];
                    v_spec_leak = [v_spec(1)*ones(Npt_cell,1); v_spec; v_spec(end)*ones(Npt_cell,1)];
                    
                    switch model_g % leaking gas properties
                        case 'simple'
                            mu_g_Lkg23  = mu_g*ones(Npt_comp,1);            % leaking ideal gas viscosity [Pa s]
                        case 'full'
                            switch model_full
                                case 'fullTP'
                                    [mu_g_Lkg23] = ThermoPhysProps('muv',T_g,p_comp,name4prop,molcomp_g,propSource,0,'PR','ChungEtal',1e-3,80); % leaking real gas properties
                                case 'fullTV'
                                    [rho_liq,rho_vap] = ThermoPhysProps('rhol,rhov',T_g,v_spec,name4prop,molcomp_g,'REFPROP_TV',0);
                                    % specific volumes of saturated liquid and vapour phase: corrective coefficient are defined inorder to avoid numerical problem with
                                    % Thermophysprops solution for saturated phase:
                                    fcorr = 0.01;
                                    v_liq = 1./rho_liq*(1-fcorr);
                                    v_vap = 1./rho_vap*(1+fcorr);
                                    % boolean vector defined to semplify the notation
                                    fact1 = quality >= 0 & quality <= 1;
                                    fact2 = quality > 1;
                                    fact3 = ~(fact1 | fact2);
                                    % leaking flow viscosity initialisation
                                    mu_g_Lkg23  = NaN(Npt_comp,1);
                                    
                                    if sum(fact1) >= 1
                                        [mu_liq] = ThermoPhysProps('mul',T_g(fact1),v_liq(fact1),name4prop,molcomp_g,'REFPROP_TV',0);
                                        [mu_vap] = ThermoPhysProps('muv',T_g(fact1),v_vap(fact1),name4prop,molcomp_g,'REFPROP_TV',0);
                                        % two-phase viscosity (EMT method by Awad et al. - 2008, see Franzetti-Persico thesys
%                                         mu_g_Lkg23(fact1) = 1/4*(sum([(3*quality(fact1)-1).*mu_vap,(3*(1-quality(fact1))-1).*mu_liq],2,'omitnan') + sqrt(sum([sum([(3*quality(fact1)-1).*mu_vap,(3*(1-quality(fact1))-1).*mu_liq],2,'omitnan').^2,8*mu_liq.*mu_vap],2,'omitnan')));
                                          mu_g_Lkg23_M1 = mu_liq.*sum([2*mu_liq, mu_vap-2.*quality(fact1).*(mu_liq-mu_vap)],2,'omitnan')./sum([2*mu_liq,mu_vap+quality(fact1).*(mu_liq-mu_vap)],2,'omitnan');          % -Maxwell 1 (gas based)
                                          mu_g_Lkg23_M2 = mu_vap.*sum([2*mu_vap,mu_liq-2.*(1-quality(fact1)).*(mu_vap-mu_liq)],2,'omitnan')./sum([2*mu_vap,mu_liq+(1-quality(fact1)).*(mu_vap-mu_liq)],2,'omitnan');  % -Maxwel 2 (liquid based)
                                          mu_g_Lkg23(fact1)    = sum([mu_g_Lkg23_M1,mu_g_Lkg23_M2],2,'omitnan')/2;
                                    end
                                    if sum(fact2) >= 1
                                        [mu_g_Lkg23(fact2)] = ThermoPhysProps('muv',T_g(fact2),v_spec(fact2),name4prop,molcomp_g,'REFPROP_TV',0);
                                    end
                                    if sum(fact3) >= 1
                                        [mu_g_Lkg23(fact3)] = ThermoPhysProps('mul',T_g(fact3),v_spec(fact3),name4prop,molcomp_g,'REFPROP_TV',0);
                                        warning('S3_Thermodynamics:viscosity','The flow is fully condensed.');
                                        SX_Logfile ('w',{strrep(lastwarn,'%','%%')});
                                    end
                            end
                    end
                    
                    [rho_o_Lkg23,~,mu_o_Lkg23,~,~] = S3_OilProperties(mean(T_lmean,'omitnan'),T_ref,rho_ref,nu_40,nu_100); % leaking oil properties
                   
                    [Fleak_g_VS_dis,Fleak_o_VS_dis,Fleak_g_VS_suc,Fleak_o_VS_suc,Fleak_g_VS_in_TV,Fleak_g_VS_in_LV,Fleak_g_VS_out,Fleak_o_VS_in_TV,Fleak_o_VS_in_LV,Fleak_o_VS_out,Fleak_o_VS_net,Dmg_VS_comp,Fleak_g_VS_cell_closed,Flow_VS] = S3_LeakageVS(MUmodel_VS,VSmodel,process,d,BI,s,VSclr,dt_vano,pos_SucClose,pos_DisOpen,V_cell,V_g,m_g,m_l,rpm,p_suc,p_comp_C,p_del,T_g,T_mix,R_g,mu_g_Lkg23,c_v,Z,rho_o_Lkg23,mu_o_Lkg23,Npt_cell,Npt_comp,dt);

                    [Fleak_g_PE_dis,Fleak_o_PE_dis,Fleak_g_PE_suc,Fleak_o_PE_suc,Fleak_g_PE_cell_closed,Fleak_o_PE_net,Fleak_g_PE_net,Dmg_PE_comp,Flow_PE] = S3_LeakagePE(MUmodel_PE,PEmodel,process,d,d_hub,PEclr,dt_vano,pos_SucOpen,pos_SucClose,pos_DisOpen,pos_DisClose,Gamma,theta_vane,V_cell,V_g,m_g,m_l,rpm,p_suc,p_comp_C,p_del,T_g,T_mix,R_g,mu_g_Lkg23,c_v,Z,rho_o_Lkg23,mu_o_Lkg23,Npt_cell,Npt_comp,dt,s);

                    if sum (NOZZLES.f_nz)>=1 && process ==1
                        Vleak_o_suc_C     = (Fleak_o_VS_suc+Fleak_o_PE_suc)/rho_o_Lkg23;   % leaking oil volume in suction [m^3/s]
                        Vleak_o_VS_out    = Fleak_o_VS_out/rho_o_Lkg23;                    % leaking oil FROM the cell [m^3/s]
                        Vleak_o_cls0_C    = (Fleak_o_VS_in_TV+Fleak_o_VS_in_LV+Fleak_o_PE_net)/rho_o_Lkg23;    % leaking oil TO the cell accounting only for oil class 0 [m^3/s]
                        Vleak_o_net       = (Fleak_o_VS_net+Fleak_o_PE_net)/rho_o_Lkg23;   % net leaking oil volume in compression [m^3/s]
                        Fleak_o_net       = Fleak_o_VS_net + Fleak_o_PE_net;
                        Fleak_o_cls0_C    = Fleak_o_VS_in_TV + Fleak_o_VS_in_LV + Fleak_o_PE_net;
                        % Temperature vector for oil leaking enthalpy flow along the compression
                        T_leak_o = [T_lmean(1)*ones(Npt_cell,1); T_lmean; T_lmean(end)*ones(Npt_cell,1)];
                        
                    elseif sum (NOZZLES.f_nz)>=1 && process ==2
                        
                        Vleak_o_suc_E     = (Fleak_o_VS_suc-Fleak_o_PE_suc)/rho_o_Lkg23;   % leaking oil volume in suction [m^3/s]
                        Vleak_o_VS_out    = Fleak_o_VS_out/rho_o_Lkg23;
                        
                        % Hy) PE oil leaking mass flow rate contriution before the first nozzle is 0
                        if nzls_suc == 0 && ~isempty(pos_nz)
                            T_lmean(1:pos_nz(1))        = 0;
                            Fleak_o_PE_net(1:pos_nz(1)) = 0;
                        end
                        
                        % Temperature vector for oil leaking enthalpy flow along the compression
                        T_leak_o = [T_lmean(1)*ones(Npt_cell,1);T_lmean; T_lmean(end)*ones(Npt_cell,1)];
                        
                        Vleak_o_cls0_E    = (Fleak_o_VS_in_TV+Fleak_o_VS_in_LV+Fleak_o_PE_net)/rho_o_Lkg23;    % leaking oil TO the cell accounting only for oil class 0 [m^3/s]
                        Vleak_o_net       = (Fleak_o_VS_net+Fleak_o_PE_net)/rho_o_Lkg23;   % net leaking oil volume in compression [m^3/s]
                        Fleak_o_net       = Fleak_o_VS_net + Fleak_o_PE_net;
                        Fleak_o_cls0_E    = Fleak_o_VS_in_TV + Fleak_o_VS_in_LV + Fleak_o_PE_net;
                        
                    else
                        rho_o_Lkg23 = ones(1,1);                % Leaking oil density [kg/m^3] (set equal to one to avoid problem in bg_4 solving the thermodynamic system)
                    end
                    
                    Niter_In = Niter_In+1;                      % Iteration counter
                    
                    % If the leakages do not converge after IterMax iteration, SVEC stop the calculation, keeping the latest result
                    if Niter_In == IterMax
                        Loop_In = 0;
                        warning( 'S3_Thermodynamics:InIter','After %d iterations, vane end-plate leakage and rotor end-plate do not converge. Result of last iteration is mantained.', IterMax);
                        SX_Logfile ('w',{lastwarn});
                    end
                end
            else % leakages disabled: no gas leakage
                Loop_In = 0; % exit outer cycle
            end  % end inner loop (vane end-plates clearance leakage)
        end
        
        clear h Al NS A B B1 Q_l C_l val pos fcorr
        %% ROTOR-STATOR LEAKAGE COMPUTATIONS %%
        if fLKG
            % definition of delivery pressure
            if isnan(p_del)
                p_out=p_comp(end);
            else
                p_out=p_del;
            end
            
            % upstream & downstream conditions
            if process == 1           % condition for compressors
                Pup   = p_out;        % upstream pressure [pa]
                Pdo   = p_suc;        % downstream pressure [pa]
                Tup_g = T_g(end);     % upstream gas temperature [K]
                Tup_o = T_lmean(end); % upstream oil temperature [K]
                v_spec_up  = v_spec(end);  % upstream specific volume [m^3/kg]
                quality_up = quality(end); % upstream quality [-]
                
            elseif process == 2       % condition for expanderers
                Pup   = p_suc;        % upstream pressure [pa]
                Pdo   = p_out;        % downstream pressure [pa]
                Tup_g = T_g(1);       % upstream gas temperature [K]
                Tup_o = T_lmean(end); % upstream oil temperature [K]
                v_spec_up  = v_spec(1);  % upstream specific volume [m^3/kg]
                quality_up = quality(1); % upstream quality [-]
            end
            
            % Leaking fluids density and viscosity
            [rho_o_RS,~,mu_o_RS,~,~] = S3_OilProperties(Tup_o,T_ref,rho_ref,nu_40,nu_100); % leaking oil properties
            switch model_g % leaking gas properties
                case 'simple'
                    rho_g_RS = Pup/(R_g*Tup_g); % leaking ideal gas density [kg/m3]
                    mu_g_RS  = mu_g;            % leaking ideal gas viscosity [Pa s]
                case 'full'
                    switch model_full
                        case 'fullTP'
                            [rho_g_RS,mu_g_RS] = ThermoPhysProps('rho,muv',Tup_g,Pup,name4prop,molcomp_g,propSource,0,'PR','ChungEtal',1e-3,80); % leaking real gas properties
                        case 'fullTV'
                            rho_g_RS = 1/v_spec_up;
                            [rho_liq_up,rho_vap_up] = ThermoPhysProps('rhol,rhov',Tup_g,v_spec_up,name4prop,molcomp_g,'REFPROP_TV',0);
                            % specific volumes of saturated liquid and vapour phase: corrective coefficient are defined inorder to avoid numerical problem with
                            % Thermophysprops solution for saturated phase:
                            fcorr = 0.01;
                            v_liq_up = 1/rho_liq_up*(1-fcorr);
                            v_vap_up = 1/rho_vap_up*(1+fcorr);
                            
                            if quality_up > 0 && quality_up < 1
                                [mu_liq_up] = ThermoPhysProps('mul',Tup_g,v_liq_up,name4prop,molcomp_g,'REFPROP_TV',0);
                                [mu_vap_up] = ThermoPhysProps('muv',Tup_g,v_vap_up,name4prop,molcomp_g,'REFPROP_TV',0);
                                % two-phase conductivity (EMT method by Awad et al. - 2008, see Franzetti-Persico thesys
                                mu_g_RS = 1/4*(sum([(3*quality_up-1).*mu_vap_up,(3*(1-quality_up)-1).*mu_liq_up],2,'omitnan') + sqrt(sum([sum([(3*quality_up-1).*mu_vap_up,(3*(1-quality_up)-1).*mu_liq_up],2,'omitnan').^2,8*mu_liq_up.*mu_vap_up],2,'omitnan')));
                            elseif quality_up >= 1
                                [mu_g_RS] = ThermoPhysProps('muv',Tup_g,v_vap_up,name4prop,molcomp_g,'REFPROP_TV',0);
                            else
                                [mu_g_RS] = ThermoPhysProps('mul',Tup_g,v_liq_up,name4prop,molcomp_g,'REFPROP_TV',0);
                                warning('S3_Thermodynamics:viscosity','The rotor-stator leaking flow is fully condensed.');
                                SX_Logfile ('w',{strrep(lastwarn,'%','%%')});
                            end
                    end
            end
            
            % Leaking flow rate
            [Fleak_g,Fleak_o,fOK] = S3_LeakageRS (MUmodel_RS,RSmodel,c,process,D,d,L,RSclr,Gamma,TgAngle,dt_vano,theta,pos_SucOpen,pos_SucClose,pos_DisOpen,pos_DisClose,pos_End,V_cell,rpm,Pup,Pdo,T_mix,R_g,rho_g_RS,mu_g_RS,n_van,V_inj_nzl,cnd_suc,c_v,rho_o_RS,mu_o_RS,Npt_cell);
                  
            %these are the variables needed to obtain the plot of leakages
            if fLKG
            Flow_RS=[Fleak_g.*ones(Npt_cell,1)./(p_suc./(R_g.*T_suc.*Z_suc)).*1e3.*60; zeros(Npt_cell*5,1); -Fleak_g.*ones(Npt_cell+1,1)./(p_suc./(R_g.*T_suc.*Z_suc)).*1e3.*60];
            if fLKG_in
            Flow_VS=[zeros(1,pos_SucClose-2*Npt_cell) ,Flow_VS./(p_suc./(R_g.*T_suc.*Z_suc)).*1e3.*60, zeros(1,length(theta_vane)-(pos_DisOpen-Npt_cell)-1) ] ;
            Flow_PE=Flow_PE./(p_suc./(R_g.*T_suc.*Z_suc)).*1e3.*60;
            end
            end
 
            % COMPRESSORS subroutine
            if sum (NOZZLES.f_nz)>=1 && process ==1       
                % in compressors the initial guess on the oil leakage flow
                % rate must be checked and possibly updated
                Vleak_o_new = Fleak_o/rho_o(1);        % mass of leaking oil via rotor-stator clearance, new iteration [m3/s]
                VleakRS_C   = NOZZLES.specs{1}.V_inj;  % mass of leaking oil via rotor-stator clearance, current iteration [m3/s]
                Tl_leak     = Tl_in(1);                % leaking oil temperature, current iteration [K]
                
                % Compute error between results of current iteration respect to previous one
                errVleak = abs((VleakRS_C-Vleak_o_new)/VleakRS_C); % leaking oil volume - error on outer iteration (RS leakage) [-]
                errTleak = abs((Tl_leak-Tup_o)/Tl_leak);           % error on leaking oil density - outer iteration (RS leakage) [-]
                
                % Check iteration error, if it is low, exit outer cycle otherwise update leaking oil mass and density and re-iterate
                if errVleak <= toll_d && errTleak <= toll_d
                    Loop_Out = 0;
                else
                    NOZZLES.specs{1}.V_inj = UndRlx*Vleak_o_new+(1-UndRlx)*VleakRS_C; % update leaking oil flow rate [m3/s]
                    Tl_in(1)               = Tup_o;                                   % update leaking oil temperature [K]
                    Niter_Out              = Niter_Out+1;                             % iteration counter
                    
                    % If the leakages do not converge after IterMax iteration,SVEC stop the calculation, keeping the latest result
                    if Niter_Out == IterMax
                        Loop_Out = 0;
                        warning( 'S3_Thermodynamics:OutIter','After %d iterations, rotor-stator leakage do not converge. Result of last iteration is mantained.', IterMax);
                        SX_Logfile ('w',{lastwarn});
                    end
                end
            
            % EXPANDERS subroutine
            elseif sum (NOZZLES.f_nz)>=1 && process == 2  
                % In expanders, the initial guess on injected oil must be
                % checked and possibly updated
                Vleak_o_new = Fleak_o/rho_o_RS; % oil leaking flow rate, new iteration [m3/s]
                
                % Compute error between results of current iteration respect to previous one
                errVleak = abs((VleakRS_E-Vleak_o_new)/VleakRS_E); % leaking oil volume - error on outer iteration (RS leakage) [-]
                
                % Check iteration error, if it is low, exit outer cycle otherwise update leaking oil mass and density and re-iterate
                if errVleak <= toll_d
                    Loop_Out = 0;
                elseif Vleak_o_new == 0 
                    Loop_Out = 0;
                    warning( 'S3_Thermodynamics:OutIter','Vleak_o_RS = 0 m3/s since no oil is present in suction.');
                    SX_Logfile ('w',{lastwarn});
                else
                    VleakRS_E = UndRlx*Vleak_o_new+(1-UndRlx)*VleakRS_E; % update leaking oil flow rate [m3/s]
                    Niter_Out = Niter_Out+1;                             % iteration counter
                    % If the leakages do not converge after IterMax iteration,SVEC stop the calculation, keeping the latest result
                    if Niter_Out == IterMax
                        Loop_Out = 0;
                        warning( 'S3_Thermodynamics:OutIter','After %d iterations, rotor-stator leakage do not converge. Result of last iteration is mantained.', IterMax);
                        SX_Logfile ('w',{lastwarn});
                    end
                end
                
            else % no oil injection, no need for iteration or correction. Only gas leakage result is needed
                Loop_Out = 0;
            end
            
        else % leakages disabled
            Fleak_g  = 0; % no oil leakage
            Fleak_o  = 0; % no gas leakage
            Loop_Out = 0; % exit outer cycle
        end  % end outer loop (rotor stator clearance leakage)
        
    end
    
    %% CALCULATION OF POWER AND WORK %%
    % Injection temperatures for each class [K]
    Tl_in_cls = arrayfun(@(x) T_l(find(~isnan(T_l(:,x)),1),x), 1:K);

    % Power evaluation
    Pind_ist(1:end-1) = -p_comp(1:end-1).*diff(V_g)/dt; % Indicated Power [W] - computed from PV
    dQprocess = sum(dQ,2,'omitnan');                    % Heat flow rate [W] for each discretization point

    % compression work evaluation
    DUgas = sum(dUgas,'omitnan');
    DUliq = sum(dUliq,'omitnan');
    DU    = DUgas + DUliq;                                  % mixture compression work - analytical [J]
    DH_g  = sum(dH_g,'omitnan');
    DH_o  = sum(dH_o,'omitnan');
    Lvano = sum(L_g,'omitnan');                            % mixture compression work - numerical [J]
    % Lvano = dt*trapz(Pind_ist(1:end-1));
    
    % compression heat evaluation
    DQ    = dt*trapz(dQprocess);                          % Heat flow rate between gas and liquid - numerical [J]
    Q     = sum(m_inj_cls.*c_l.*(T_l(end,:)-Tl_in_cls));  % Heat flow rate between gas and liquid - analytical [J]
    Q_th  = Q*c*n_van*rpm/60;                             % Thermal power [W]
    clear Tl_in_cls dQprocess DUgas DUliq dQ

    %% COMPARISON WITH THE POLYTROPIC PROCESS %%
    cx(1:end-1)     = c_v-dt*Pind_ist(1:end-1)./(diff(T_g))./m_g(1:end-1);   % Specific heat for the polytropic process [J/kgK]
    n_pltr(1:end-1) = (cx(1:end-1)-(c_v+R_g))./(cx(1:end-1)-c_v);            % Polytropic index [adim]

    %% REMOVING FICTIOUS LEAKAGE NOZZLE %%
    % only output variables are restored. Commented variables are affected
    % by fictious nozzle, are not more used in SVEC.
    if fLKG && process==1 && sum(NOZZLES.f_nz)>=1
        M         = M-1;               % restore original nozzle number
        V_inj_nzl = V_inj_nzl(2:end);  % restore original nozzle injected volume [m3]
        m_inj_nzl = m_inj_nzl(2:end);  % restore original nozzle injected mass [kg]
        V_inj_cls = V_inj_cls(2:end);  % restore original class injected mass [kg]
        rho_o     = rho_o(2:end);      % restore original oil class density [kg/m3]
        nu_o      = nu_o(2:end);       % restore original oil viscosity [cSt]
        sigma_o   = sigma_o(2:end);    % restore original oil surface tension [N/m]
        Dg        = Dg(2:end);         % restore original class droplet diameter [m]
        fg        = fg(2:end);         % restore original class frequency [-]
        T_l       = T_l(:,2:end);      % restore original class temperature [-]
        % K       = K-1;               % restore original number of classes
        % mu_o    = nu_o(2:end);       % restore original oil viscosity [Pa s]
    end
    

    %% FLOW RATES %%
    if process == 1
        s_lkg = -1;
    else
        s_lkg = 1;
    end
    
    m_gas_Id = p_suc*V_comp(1)/(R_g*T_suc*Z_suc)*c*n_van*rpm/60;                                      % ideal gas mass flow rate [kg/s]
    m_gas    = m_g(end)*c*n_van*rpm/60 + s_lkg*Fleak_g + s_lkg*Fleak_g_VS_dis + Fleak_g_PE_dis; % real gas mass flow rate [kg/s]
    q_gas_Id = (m_gas_Id)/(p_suc/(R_g*T_suc*Z_suc)); 
    

    % ideal gas volumetric flow rate [m3/s]
    q_gas    = (m_gas)/(p_suc/(R_g*T_suc*Z_suc));                                               % real gas volumetric flow rate [m3/s]
    m_liq_Id = sum(m_inj_nzl)*c*n_van*rpm/60;                                                   % mass flow rate of injected oil [kg/s]
    m_liq    = sum(m_cell_cls(end,1:end),'omitnan') + s_lkg*Fleak_o + s_lkg*Fleak_o_VS_dis + Fleak_o_PE_dis;                 % real oil mass flow rate [kg/s]
    epsm     = m_inj_nzl./mean(m_g);
   
    % epsm for each injector
    % =====================================================================
    % =====================================================================
    Q_g_RS_suction = Fleak_g/(p_suc/(R_g*T_suc*Z_suc));
    Q_g_PE_suction = Fleak_g_PE_suc/(p_suc/(R_g*T_suc*Z_suc));
    Q_g_VS_suction = Fleak_g_VS_suc/(p_suc/(R_g*T_suc*Z_suc));
    Q_g_PE_closed_cell = Fleak_g_PE_cell_closed/(p_suc/(R_g*T_suc*Z_suc));
    Q_g_VS_closed_cell = Fleak_g_VS_cell_closed/(p_suc/(R_g*T_suc*Z_suc));
    Q_g_RS_dis = Fleak_g/(p_suc/(R_g*T_suc*Z_suc));
    Q_g_PE_dis = Fleak_g_PE_dis/(p_suc/(R_g*T_suc*Z_suc));
    Q_g_VS_dis = Fleak_g_VS_dis/(p_suc/(R_g*T_suc*Z_suc));
    % =====================================================================
    % =====================================================================    
      
    %% CHECKS %%
    Bi_max  = max(Bi);                                          % Check on Biot number
    errDU   = abs(DU-(sum([Lvano,DH_g,DH_o],'omitnan')))/DU;    % Verify the first law of thermodynamics DU = Q + L with Q=0 (adiabatic system)
    m_in    = p_comp(1)*V_g(1)/R_g/T_g(1)/Z(1) + s_lkg*Fleak_g_VS_suc*dt_vano - Fleak_g_PE_suc*dt_angle_suc;       % Verify the mass balance m_in = m_out (steady state)
    m_out   = p_comp(end)*V_g(end)/R_g/T_g(end)/Z(end) + s_lkg*Fleak_g_VS_dis*dt_vano + Fleak_g_PE_dis*dt_angle_dis;
     
    
    errmass = abs(m_in - m_out)/m_in;            % error on mass balance
    errQ    = abs(Q-DQ)/DQ;                      % Check on heat transfer

    if sum(m_g <0) >= 1
        warning( 'S3_Thermodynamics:GasMass','Gas mass in the cell is negative. Try a more conservative leakage model');
        SX_Logfile ('w',{lastwarn});
    end
    
    if sum(m_l <0) >= 1
        warning( 'S3_Thermodynamics:OilMass','Oil mass in the cell is negative. Try a more conservative leakage model');
        SX_Logfile ('w',{lastwarn});
    end
    
    if m_gas < 0
        warning( 'S3_Thermodynamics:GasFlow','Gas flow rate is negative. Try a more conservative leakage model');
        SX_Logfile ('w',{lastwarn});
    end
    
    if m_liq < 0
        warning( 'S3_Thermodynamics:OilFlow','Oil flow rate is negative. Try a more conservative leakage model');
        SX_Logfile ('w',{lastwarn});
    end
    
    if Bi_max > 0.1
        warning( 'S3_Thermodynamics:numeroBiot','Model cannot be considered with concentrated parameter. Max Biot number is: %2.2d',max(Bi));
        SX_Logfile ('w',{lastwarn});
    end

    if abs(errDU) > toll_d
        warning('S3_Thermodynamics:energia','System: input and output energy do no match. Estimated relative error: %2.4f %%', abs(errDU)*100);
        SX_Logfile ('w',{strrep(lastwarn,'%','%%')});
    end

    if abs(errmass) > toll_d
        warning('S3_Thermodynamics:massa','System: input and output mass do not match. Estimated relative error : %2.4f %%', abs(errmass)*100);
        SX_Logfile ('w',{strrep(lastwarn,'%','%%')});
    end

    if abs(errQ) > toll_d
        warning('S3_Thermodynamics:calore','System: oil exchanged heat not matching. Estimated relative error: %2.4f %%', abs(errQ)*100);
        SX_Logfile ('w',{strrep(lastwarn,'%','%%')});
    end
    clear Bi_max errDU errmass errQ m_in m_out Z
end


function F=system(T,Ki,T_g,T_l,Tleak_g_23,Tleak_o_23,T_0,p_0,v_spec_0,p_comp,p_leak_23,rho_o_Lkg23,Vg,m_g,v_spec,v_spec_leak_23,m_cell_cls,Fleak_g_23,Fleak_o_23,c_l,h,Al,NS,dt,model_full,name4prop,z,source)
% Call ThermoPhysProps for points (i) and (i+1)
switch model_full
    case 'fullTP'
        [u_g]    = ThermoPhysProps('u',[T_g,T(1)],[p_comp,p_comp],name4prop,z,source,0,'PR','ChungEtal',1e-3,80);
        [h_g]    = ThermoPhysProps('h',[Tleak_g_23,T_0],[p_leak_23,p_0],name4prop,z,source,0,'PR','ChungEtal',1e-3,80);
    case 'fullTV'
        [u_g]    = ThermoPhysProps('u',[T_g,T(1)],[v_spec(1),v_spec(2)],name4prop,z,'REFPROP_TV',0);
        [h_g]    = ThermoPhysProps('h',[Tleak_g_23,T_0],[v_spec_leak_23,v_spec_0],name4prop,z,'REFPROP_TV',0);
end

% Build equation system
F        = NaN(Ki+1,1);
DH_g_n   = (Fleak_g_23(1)*(h_g(1)-h_g(5)) + Fleak_g_23(2)*(h_g(2)-h_g(5)) - Fleak_g_23(3)*(h_g(3)-h_g(5)))*dt + Fleak_g_23(4)*dt*(h_g(4)-h_g(5));
DH_o_n   =(Fleak_o_23(1)*(Tleak_o_23(1)-T_0) + Fleak_o_23(2)*(Tleak_o_23(2)-T_0) - Fleak_o_23(3)*(Tleak_o_23(3)-T_0))*c_l*dt + Fleak_o_23(4)*c_l*dt*(Tleak_o_23(4)-T_0) + ...
          + p_leak_23(1)*Fleak_o_23(1)*dt/rho_o_Lkg23 + p_leak_23(2)*Fleak_o_23(2)*dt/rho_o_Lkg23 - p_leak_23(3)*Fleak_o_23(3)*dt/rho_o_Lkg23 + p_leak_23(4)*Fleak_o_23(4)*dt/rho_o_Lkg23;
        
% System
F(1)     = ((m_g.*(u_g(2)-u_g(1))+sum(m_cell_cls.*c_l.*(T(2:end)-T_l)))./(-p_comp.*(Vg(2)-Vg(1))+DH_g_n+DH_o_n))-1;
F(2:end) = (m_cell_cls'.*c_l.*(T(2:end)'-T_l'))./(-h'.*Al'.*NS'.*(T(2:end)'-T(1)).*dt)-1;
end