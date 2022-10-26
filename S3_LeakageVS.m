function [Fleak_g_VS_dis,Fleak_o_VS_dis,Fleak_g_VS_suc,Fleak_o_VS_suc,Fleak_g_VS_in_TV,Fleak_g_VS_in_LV,Fleak_g_VS_out,Fleak_o_VS_in_TV,Fleak_o_VS_in_LV,Fleak_o_VS_out,Fleak_o_VS_net,Dmg_VS_comp,Fleak_g_VS_net,Flow_VS] = S3_LeakageVS(MUmodel_VS,VSmodel,process,d,BI,s,VSclr,dt_vano,pos_SucClose,pos_DisOpen,V_cell,V_g,m_g,m_o,rpm,p_suc,p_comp,p_del,T_g,T_mix,R_g,mu_g_VS,c_v,Z,rho_o_VS,mu_o_VS,Npt_cell,Npt_comp,dt)
% This function compute leakages through clearance between vane and end-plates (path 3)
%
% INPUT
% MUmodel_VS [string]  : viscosity model for two-phase gas-liquid mixture
% VSmodel [string]  : vane-end-plate leakage model
% process           : process selected (1: compression  2: expansion
% d [m]             : rotor diameter
% BI [m]            : real vane excursion measured on vane axis
% s [m]             : vane thickness
% VSclr [m]         : vane-end-plate clearance size
% dt_vano [s]       : time needed for a vane to sweep a cell [s]
% pos_SucClose [-]  : discretized suction close angle position in array theta
% pos_DisOpen [-]   : discretized delivery open angle position in array theta
% V_cell [m^3]      : cell volumes array of a complete vane revolution (-gamma:2pi)
% V_g    [m^3]      : gas volume for each discretization point
% m_g    [kg]       : gas mass for each discretization point
% m_o    [kg]       : oil mass for each discretization point
% rpm [rpm]         : rotor angular velocity
% p_suc [Pa]        : suction pressure
% p_comp [Pa]       : vector of pressures from suction closure to discharge aperture
% p_del  [Pa]       : delivery pressure
% T_g [K]           : temperature array
% T_mix [K]         : oil + gas adiabatic mixing temperature for each discretization step
% R_g [J/kgK]       : specific gas constant
% mu_g_VS [Pa s]    : leaked gas dynamic viscosity
% c_v [J/kg K]      : specific heat at costant volume
% Z   [-]              : compressibility factor 
% rho_o_VS [kg/m3]  : leaked oil density
% mu_o_VS [Pa s]    : leaked oil kinematic viscosity
% Npt_cell [-]      : number of discretization points between 2 vanes
% Npt_comp [-]      : discretization points of closed cell process 
% dt [s]            : time step 

% OUTPUT
% Fleak_g_VS_dis [kg/s]     : leaking gas mass flow rate in discharge through path 3
% Fleak_o_VS_dis [kg/s]     : leaking oil mass flow rate in discharge through path 3
% Fleak_g_VS_suc [kg/s]     : leaking gas mass flow rate in suction through path 3
% Fleak_o_VS_suc [kg/s]     : leaking oil mass flow rate in suction through path 3
% Fleak_g_VS_in_TV [kg/s]   : leaking gas mass flow rate TO the cell through path 3 on trailing vane (TV)
% Fleak_g_VS_in_LV [kg/s]   : leaking gas mass flow rate TO the cell through path 3 on leading vane (LV)
% Fleak_g_VS_out   [kg/s]   : leaking gas mass flow rate FROM the cell through path 3 
% Fleak_o_VS_in_TV [kg/s]   : leaking oil mass flow rate TO the cell through path 3 on trailing vane (TV)
% Fleak_o_VS_in_LV [kg/s]   : leaking oil mass flow rate TO the cell through path 3 on leading vane (LV)
% Fleak_o_VS_out   [kg/s]   : leaking oil mass flow rate FROM the cell through path 3 
% Fleak_o_VS_net [kg/s]     : net leaking oil mass flow rate along compression process through path 3 
% Dmg_VS_comp [kg]          : cumulated gas mass variation due to leakage through path 3 along compression process

%NOTE:  Franzetti-Persico theoretical description and SVEC numerical application might sound different.   
%This is only due to the different references of the theoretical and numerical approach, but results are coherent.
%For further information, please read comments on Franzetti-Persico thesis (page 72), available both on Mattei servers and in the shared Onedrive folder.   

    %% DEFINITIONS %%

    % GEOMETRY
    r = d/2;             % rotor radius  [m]
    
    % THERMODYNAMICS
    c_p   = c_v+R_g;     % gas specific heat at costant pressure [J/kgK]
    Gamma = c_p/c_v;     % specific heat ratio [-]
         
    % PRECEDING cell    
    % Pressure
    p_minus = [p_suc*ones(Npt_cell,1);p_comp(1:end-Npt_cell)];   % following cell pressure [Pa]    
    
    % FOLLOWING cell    
    % Pressure
    if isnan(p_del)                                            % definition of delivery pressure
        p_out = p_comp(end);
    else
        p_out = p_del;
    end
    p_plus = [p_comp(Npt_cell+1:end);p_out*ones(Npt_cell,1)];   % preceding cell pressure [Pa]
    
    % CURRENT cell
    % Compressibility factor
    Z_cell   = Z;
    
    % Temperature
    Tg_cell  = T_g;                     % gas temperature in the current cell [K]
    T_b_cell = T_mix;                  % leakage temperature from the current cell [K]
    
    % Mass
    m_g_cell = m_g;                    % gas mass in the current cell [kg]
    m_o_cell = m_o;                    % liquid mass in the current cell [kg]
    
    % Mass fraction
    X_g_cell = m_g_cell./(m_g_cell+m_o_cell);      % gas mass fraction in the current cell [-]
    
    % Volume fraction
    X_vol_cell = V_g./V_cell(pos_SucClose:pos_DisOpen-Npt_cell)';  % gas volume fraction in the current cell [-]
    
    % Dynamic viscosity
    mu_g_VS_cell = mu_g_VS;
    
    % Density
    rho_g_cell = p_comp./(Z_cell.*R_g.*Tg_cell);         % gas density in the current cell [kg/m^3]
    rho_b_cell = 1./sum([X_g_cell./rho_g_cell,(1-X_g_cell)./rho_o_VS],2,'omitnan');    % bulk density in the current cell [kg/m^3]
    
    % UPSTREAM and DOWNSTREAM leakage parameters definition
    switch process
        case 1                         % compression process
            % Geometry
            W_cl_cell = BI(pos_SucClose-Npt_cell:pos_DisOpen-2*Npt_cell)';                  % real trailing vane excursion measured on vane axis
            W_cl_up   = BI(pos_SucClose-Npt_cell+Npt_cell:pos_DisOpen-Npt_cell)';           % real leading vane excursion measured on vane axis
            
            % Pressures
            p_do = p_minus;            % leakage downstream pressure [Pa]
            p_up = p_plus;             % leakage upstream pressure [Pa]
            
            % Compressibility factor
            Z_up   = [Z(Npt_cell+1:end);Z(end)*ones(Npt_cell,1)];  
            
            % Temperatures
            Tg_up  = [T_g(Npt_cell+1:end);T_g(end)*ones(Npt_cell,1)];          % gas temperature in the preceding cell [K]
            T_b_up = [T_mix(Npt_cell+1:end);T_mix(end)*ones(Npt_cell,1)];      % leakage upstream temperature  [K]
            
            % Mass    
            m_g_up = [m_g(Npt_cell+1:end);m_g(end)*ones(Npt_cell,1)];          % gas mass in the preceding cell [kg]
            m_o_up = [m_o(Npt_cell+1:end);m_o(end)*ones(Npt_cell,1)];          % liquid mass in the preceding cell [kg]
            
            % Dynamic viscosity
            mu_g_VS_up = [mu_g_VS(Npt_cell+1:end);mu_g_VS(end)*ones(Npt_cell,1)];

            % Volume fraction    
            X_vol_up = [V_g(Npt_cell+1:end); V_g(end)*ones(Npt_cell,1)]./[V_cell(pos_SucClose+Npt_cell:pos_DisOpen-Npt_cell)';V_cell(pos_DisOpen-Npt_cell)*ones(Npt_cell,1)];           % gas volume fraction in the current cell [-]
                                                 
        case 2                         % expansion process
            % Geometry
            W_cl_cell = BI(pos_SucClose-Npt_cell+Npt_cell:pos_DisOpen-Npt_cell)';     % real leading vane excursion measured on vane axis
            W_cl_up   = BI(pos_SucClose-Npt_cell:pos_DisOpen-2*Npt_cell)';            % real trailing vane excursion measured on vane axis
            
            % Pressures
            p_do = p_plus;             % leakage downstream pressure [Pa]
            p_up = p_minus;            % leakage upstream pressure [Pa]
            
            % Compressibility factor
            Z_up   = [Z(1)*ones(Npt_cell,1);Z(1:end-Npt_cell)];  
            
            % Temperatures
            Tg_up  = [T_g(1)*ones(Npt_cell,1);T_g(1:end-Npt_cell)];          % gas temperature in the preceding cell [K]
            T_b_up = [T_mix(1)*ones(Npt_cell,1);T_mix(1:end-Npt_cell)];      % leakage upstream temperature  [K]
            
            % Mass    
            m_g_up = [m_g(1)*ones(Npt_cell,1);m_g(1:end-Npt_cell)];          % gas mass in the preceding cell [kg]
            m_o_up = [m_o(1)*ones(Npt_cell,1);m_o(1:end-Npt_cell)];          % liquid mass in the preceding cell [kg]
            
            % Dynamic viscosity
            mu_g_VS_up = [mu_g_VS(1)*ones(Npt_cell,1);mu_g_VS(1:end-Npt_cell)];

            % Volume fraction    
            X_vol_up = [V_g(1)*ones(Npt_cell,1);V_g(1:end-Npt_cell)]./[V_cell(pos_SucClose)*ones(Npt_cell,1);V_cell(pos_SucClose:pos_DisOpen-2*Npt_cell)'];           % gas volume fraction in the current cell [-]
    end
    
    % Mass fraction   
    X_g_up = m_g_up./(m_g_up+m_o_up);   % gas mass fraction in the preceding cell [-]
   
    % Density    
    rho_g_up = p_up./(Z_up.*R_g.*Tg_up);                                  % gas density in the preceding cell [kg/m^3]    
    rho_b_up = 1./sum([X_g_up./rho_g_up,(1-X_g_up)./rho_o_VS],2,'omitnan');            % bulk density in the preceding cell [kg/m^3]

    % Viscosity model
    switch MUmodel_VS
        case 'Dukler'
            mu_b_cell  = rho_b_cell.*sum([X_g_cell.*mu_g_VS_cell./rho_g_cell,(1-X_g_cell).*mu_o_VS./rho_o_VS],2,'omitnan');
            mu_b_up    = rho_b_up.*sum([X_g_up.*mu_g_VS_up./rho_g_up,(1-X_g_up).*mu_o_VS./rho_o_VS],2,'omitnan');            
        case 'Maxwell'
            mu_b_M1_cell = mu_o_VS.*sum([2*mu_o_VS*ones(Npt_comp,1),mu_g_VS_cell-2.*X_g_cell.*(mu_o_VS-mu_g_VS_cell)],2,'omitnan')./sum([2*mu_o_VS*ones(Npt_comp,1),mu_g_VS_cell+X_g_cell.*(mu_o_VS-mu_g_VS_cell)],2,'omitnan');          % -Maxwell 1 (gas based)
            mu_b_M2_cell = mu_g_VS_cell.*sum([2*mu_g_VS_cell,mu_o_VS-2.*(1-X_g_cell).*(mu_g_VS_cell-mu_o_VS)],2,'omitnan')./sum([2*mu_g_VS_cell,mu_o_VS+(1-X_g_cell).*(mu_g_VS_cell-mu_o_VS)],2,'omitnan');  % -Maxwel 2 (liquid based)
            mu_b_cell    = sum([mu_b_M1_cell,mu_b_M2_cell],2,'omitnan')/2;
            
            mu_b_M1_up = mu_o_VS.*sum([2*mu_o_VS*ones(Npt_comp,1),mu_g_VS_up-2.*X_g_up.*(mu_o_VS-mu_g_VS_up)],2,'omitnan')./sum([2*mu_o_VS*ones(Npt_comp,1),mu_g_VS_up+X_g_up.*(mu_o_VS-mu_g_VS_up)],2,'omitnan');            % -Maxwell 1 (gas based)
            mu_b_M2_up = mu_g_VS_up.*sum([2*mu_g_VS_up,mu_o_VS-2.*(1-X_g_up).*(mu_g_VS_up-mu_o_VS)],2,'omitnan')./sum([2*mu_g_VS_up,mu_o_VS+(1-X_g_up).*(mu_g_VS_up-mu_o_VS)],2,'omitnan');    % -Maxwel 2 (liquid based)
            mu_b_up    = sum([mu_b_M1_up,mu_b_M2_up],2,'omitnan')/2;  
            
            if sum(isnan(mu_b_M1_cell)) >= 1
                mu_b_cell =  mu_b_M2_cell;
            end
            if sum(isnan(mu_b_M1_up)) >= 1
                mu_b_up =  mu_b_M2_up;
            end 
        case 'Awad'
            mu_b_cell   = 1/4*(sum([(3*X_g_cell-1).*mu_g_VS,(3*(1-X_g_cell)-1).*mu_o_VS],2,'omitnan') + sqrt(sum([sum([(3*X_g_cell-1).*mu_g_VS,(3*(1-X_g_cell)-1).*mu_o_VS],2,'omitnan').^2,8*mu_o_VS.*mu_g_VS],2,'omitnan')));
            mu_b_up     = 1/4*(sum([(3*X_g_up-1).*mu_g_VS,(3*(1-X_g_up)-1).*mu_o_VS],2,'omitnan') + sqrt(sum([sum([(3*X_g_up-1).*mu_g_VS,(3*(1-X_g_up)-1).*mu_o_VS],2,'omitnan').^2,8*mu_o_VS.*mu_g_VS],2,'omitnan')));
    end

    %% MODEL SELECTION %%
    % INITIALIZATION
    FlowLeak_cell    = NaN(Npt_comp,1);
    FlowLeak_up      = NaN(Npt_comp,1);
    Fleak_g_VS_in_TV = zeros(Npt_comp,1);
    Fleak_g_VS_in_LV = zeros(Npt_comp,1);
    Fleak_g_VS_out   = zeros(Npt_comp,1);
    Fleak_o_VS_in_TV = zeros(Npt_comp,1);
    Fleak_o_VS_in_LV = zeros(Npt_comp,1);
    Fleak_o_VS_out   = zeros(Npt_comp,1);

    switch VSmodel
        % FlowLeak_VS : [kg/s] net leaked mass flow rate through path 3 for each discretization point along compression (>0 entering; <0 exiting the cell)
        case 'Yuan'
            [FlowLeak_cell] = Yuan(p_comp,p_do,T_b_cell,R_g,mu_b_cell,rho_b_cell,Gamma,Z_cell,s,VSclr,W_cl_cell); % Current cell
            [FlowLeak_up]   = Yuan(p_up,p_comp,T_b_up,R_g,mu_b_up,rho_b_up,Gamma,Z_up,s,VSclr,W_cl_up);         % Preceding cell
        case 'Poiselle-Couette'
            [FlowLeak_cell] = PoiselleCouette(p_comp,p_do,rho_b_cell,mu_b_cell,process,rpm,s,VSclr,W_cl_cell,r); % Current cell
            [FlowLeak_up]   = PoiselleCouette(p_up,p_comp,rho_b_up,mu_b_up,process,rpm,s,VSclr,W_cl_up,r);       % Preceding cell
        case 'Badr'
            [FlowLeak_cell] = Badr (p_comp,p_do,T_b_cell,rho_b_cell,c_v,R_g,VSclr,W_cl_cell); % Current cell
            [FlowLeak_up]   = Badr (p_up,p_comp,T_b_up,rho_b_up,c_v,R_g,VSclr,W_cl_up); % Preceding cell
        case 'Ishii'
            [FlowLeak_cell] = Ishii (p_comp,p_do,rho_b_cell,mu_b_cell,s,VSclr,W_cl_cell); % Current cell
            [FlowLeak_up]   = Ishii (p_up,p_comp,rho_b_up,mu_b_up,s,VSclr,W_cl_up); % Preceding cell
        case 'Suefuji'
            [FlowLeak_cell] = Suefuji(p_comp,p_do,T_b_cell,VSclr,W_cl_cell,mu_b_cell,R_g,Gamma); % Current cell
            [FlowLeak_up]   = Suefuji(p_up,p_comp,T_b_up,VSclr,W_cl_up,mu_b_up,R_g,Gamma); % Preceding cell
        case 'Yanagisawa'
            [FlowLeak_cell] = Yanagisawa(p_comp,p_do,T_b_cell,c_v,s,VSclr,W_cl_cell,mu_b_cell,R_g); % Current cell
            [FlowLeak_up]   = Yanagisawa(p_up,p_comp,T_b_up,c_v,s,VSclr,W_cl_up,mu_b_up,R_g); % Preceding cell
    end
 
    
    % Leakage contribution from gas and oil through path 3 for each discretization point along compression [kg/s]
    FlowLeak_cell_gas = FlowLeak_cell.*X_vol_cell;
    FlowLeak_cell_oil = FlowLeak_cell.*(1-X_vol_cell);

    FlowLeak_up_gas =  FlowLeak_up.*X_vol_up;
    FlowLeak_up_oil = FlowLeak_up.*(1-X_vol_up);
    
    % Terms that often appear in formulas
    fact1 = FlowLeak_cell_gas > 0;
    fact2 = FlowLeak_up_gas > 0;
    fact3 = FlowLeak_cell_oil > 0;
    fact4 = FlowLeak_up_oil > 0;
    
    switch process
        case 1          % compression process
        % SUCTION
        FlowLeak_Suc_gas = FlowLeak_cell_gas(1:Npt_cell);   % Gas leaking flow rate in suction during open cell phase [kg/s] (SEE NOTE)
        FlowLeak_Suc_oil = FlowLeak_cell_oil(1:Npt_cell);   % Oil leaking flow rate in suction during open cell phase [kg/s] (SEE NOTE)
        
        % DISCHARGE
        FlowLeak_Dis_gas = FlowLeak_up_gas(end-Npt_cell+1:end);  % Gas leaking flow rate in discharge during open cell phase [kg/s] (SEE NOTE)
        FlowLeak_Dis_oil = FlowLeak_up_oil(end-Npt_cell+1:end);  % Oil leaking flow rate in discharge during open cell phase [kg/s] (SEE NOTE)
       
       
        % IN-OUT leakages in current cell
        Fleak_g_VS_in_TV(~fact1)  = - FlowLeak_cell_gas(~fact1);
        Fleak_g_VS_in_LV(fact2)   = FlowLeak_up_gas(fact2);
        Fleak_g_VS_out(fact1)     = FlowLeak_cell_gas(fact1);
        Fleak_g_VS_out(~fact2)    = Fleak_g_VS_out(~fact2) - FlowLeak_up_gas(~fact2);
        
        Fleak_o_VS_in_TV(~fact3)  = - FlowLeak_cell_oil(~fact3);
        Fleak_o_VS_in_LV(fact4)   = FlowLeak_up_oil(fact4);
        Fleak_o_VS_out(fact3)     = FlowLeak_cell_oil(fact3);
        Fleak_o_VS_out(~fact4)    = Fleak_o_VS_out(~fact4) - FlowLeak_up_oil(~fact4);
        
         %COMPRESSION PHASE
        Flow_VS=[FlowLeak_Suc_gas',(Fleak_g_VS_in_TV+Fleak_g_VS_in_LV-Fleak_g_VS_out)',  -FlowLeak_Dis_gas'];
   
        case 2          % expansion process
        % SUCTION
        FlowLeak_Suc_gas = FlowLeak_up_gas(1:Npt_cell);   % Gas leaking flow rate in suction during open cell phase [kg/s] (SEE NOTE)
        FlowLeak_Suc_oil = FlowLeak_up_oil(1:Npt_cell);   % Oil leaking flow rate in suction during open cell phase [kg/s] (SEE NOTE)

        % DISCHARGE
        FlowLeak_Dis_gas = FlowLeak_cell_gas(end-Npt_cell+1:end);  % Gas leaking flow rate in discharge during open cell phase [kg/s] (SEE NOTE)
        FlowLeak_Dis_oil = FlowLeak_cell_oil(end-Npt_cell+1:end);  % Oil leaking flow rate in discharge during open cell phase [kg/s] (SEE NOTE)
        
       
        % IN-OUT leakages in current cell
        Fleak_g_VS_in_TV(fact2)    = FlowLeak_up_gas(fact2);
        Fleak_g_VS_in_LV(~fact1)   = - FlowLeak_cell_gas(~fact1);
        Fleak_g_VS_out(fact1)      = FlowLeak_cell_gas(fact1);
        Fleak_g_VS_out(~fact2)     = Fleak_g_VS_out(~fact2) - FlowLeak_up_gas(~fact2);
        
        Fleak_o_VS_in_TV(fact4)    = FlowLeak_up_oil(fact4);
        Fleak_o_VS_in_LV(~fact3)   = - FlowLeak_cell_oil(~fact3);
        Fleak_o_VS_out(fact3)      = FlowLeak_cell_oil(fact3);
        Fleak_o_VS_out(~fact4)     = Fleak_o_VS_out(~fact4) - FlowLeak_up_oil(~fact4);
        
        %EXPANSION PHASE
        Flow_VS=[-FlowLeak_Suc_gas',(Fleak_g_VS_in_TV+Fleak_g_VS_in_LV-Fleak_g_VS_out)',  FlowLeak_Dis_gas'];
   
    end

    %% COMPUTATION

    Fleak_g_VS_dis = sum(FlowLeak_Dis_gas*dt)/(dt_vano);  % leaking gas mass flow rate in discharge through path 3 [kg/s]
    Fleak_o_VS_dis = sum(FlowLeak_Dis_oil*dt)/(dt_vano);  % leaking oil mass flow rate in discharge through path 3 [kg/s]
    Fleak_g_VS_suc = sum(FlowLeak_Suc_gas*dt)/(dt_vano);  % leaking gas mass flow rate in suction through path 3 [kg/s]
    Fleak_o_VS_suc = sum(FlowLeak_Suc_oil*dt)/(dt_vano);  % leaking oil mass flow rate in suction through path 3

    % Fleak_g_VS_net = FlowLeak_up_gas-FlowLeak_cell_gas;  % net leaking gas mass flow rate along compression process through path 3 [kg/s]
    Fleak_o_VS_net = FlowLeak_up_oil-FlowLeak_cell_oil;  % net leaking oil mass flow rate along compression process through path 3 [kg/s]

    % Dmg_VS_comp = cumsum(Fleak_g_VS_net*dt); % cumulated gas mass variation due to leakage through path 3 along compression process [kg]
    % Dmo_VS_comp = cumsum(Fleak_o_VS_net*dt); cumulated oil mass variation due to leakage through path 3 along compression process [kg]
    Dmg_VS_comp = cumsum((Fleak_g_VS_in_TV+Fleak_g_VS_in_LV-Fleak_g_VS_out)*dt);
        
   
    Fleak_g_VS_net=Dmg_VS_comp(end)/(dt_vano);
    
    % removing the last element since correspond to the opening of the delivery
    Fleak_o_VS_net(end)   = [];
    Fleak_g_VS_in_TV(end) = [];
    Fleak_g_VS_in_LV(end) = [];
    Fleak_g_VS_out(end)   = [];
    Fleak_o_VS_in_TV(end) = [];
    Fleak_o_VS_in_LV(end) = [];
    Fleak_o_VS_out(end)   = [];                  

end


function [Mleak] = Yuan (Pup,Pdo,Tup,R_g,mu_b,rho_b,Gamma,Z_g,L_cl,h_cl,W_cl)
%This function compute leakages through axial clearance between vane side
%and end plates (path 3), using Yuan approach

% INPUT
% Pup [pa]       : upstream pressure
% Pdo [pa]       : downstream pressure
% T_up [K]       : upstream bulk temperature   
% r [m]          : rotor diameter
% R_g [J/kgK]    : specific gas constant
% mu_b [Pa s]    : bulk dynamic viscosity
% rho_b [kg/m^3] : bulk density
% Gamma [-]      : specific heat ratio
% L_cl [m]       : clearance length
% h_cl [m]       : clearance height
% W_cl [m]       : clearance width
%
% OUTPUT
% Mleak [kg/s] : leakage mass flow rate through path 3
%
% NOTES: - subscripts stand for:
%              up : upstream of clearance
%              do : downstream of clearance
% HP: - Incompresible flow
%     - Real gas
%     - Adiabatic Flow
%     - Viscous and Inertia forces are considered at the same time
%     - Constant clearance heigth
%     - Linear pressure dinstribution along the path
%     
% Model inspired by: "CALCULATING MODEL AND EXPERIMENTAL INVESTIGATION OF GAS LEAKAGE ", X. Yuan et Al., 1992
% See Franzetti-Persico 2019 thesis for more info

    %% COMPUTATION %%
    a = (Pdo-Pup)./L_cl.*6.*Z_g./(5.*Gamma.*R_g.*Tup.*h_cl.^2);
    b = -12*mu_b.*Z_g./h_cl.^3;
    c = -(Pdo-Pup)./L_cl;
    Delta = b.^2-4.*a.*c;
    x = (-b-sqrt(Delta))./(2*a);  %considering only positive solution
    x(isnan(x)) = 0;

    Mleak = rho_b.*W_cl.*x; % leakage mass flow rate through path 3 [kg/s]   
end

function [Mleak] = PoiselleCouette (Pup,Pdo,rho_b,mu_b,process,rpm,L_cl,h_cl,W_cl,r)
%This function compute leakages through axial clearance between vane side
%and end plates (path 3), using poiselle-Couette approach for 2D flow
% 
% INPUT
% Pup [pa]     : upstream pressure
% Pdo [pa]     : downstream pressure
% rho_b        : upstream bulk density
% mu_b [Pa s]  : bulk dynamic viscosity
% process      : process selected (1: compression  2: expansion
% rpm [rpm]    : rotor velocity
% L_cl [m]       : clearance length
% h_cl [m]       : clearance height
% W_cl [m]       : clearance width
% r [m]        : rotor diameter

% OUTPUT
% Mleak [kg/s] : leakage mass flow rate through path 3
%  
% HP: - Uncompresible flow
%     - Laminar flow
%     - Viscous flow
%     - Poiselle-Couette flow (wall velocity accounted)
%
% Model inspired by: "Simulation and comparison of leakage characteristics of R290 in rolling piston type rotary compressor", Dehua Cai, 2015
% See Franzetti-Persico 2019 thesis for more info

    %% DEFINITIONS %%
    Uwall = rpm/60*2*pi; % rotor surface velocity [rad/s]

    if process == 1 % compressione
        s = +1;
    else            % espansione
        s = -1;
    end

    %% CALCULATIONS %%
    Mleak_P = 1./(12.*mu_b.*L_cl./h_cl.^3).*W_cl.*(Pup-Pdo).*rho_b;  % leakage mass flow rate due to pressure gradient [kg/s]
    Mleak_V = W_cl.*rho_b.*h_cl./2*r.*Uwall;                         % leakage mass flow rate due to wall velocity [kg/s]
    Mleak   = Mleak_P + s * Mleak_V;                                 % total leakage mass flow rate through path 3 [kg/s]
end

function [Mleak] = Badr (Pup,Pdo,T_b,rho_b,c_v,R_g,h_cl,W_cl)
%This function compute leakages through axial clearance between vane side
%and end plates (path 3), using Badr compressible approach
% 
% INPUT
% Pup [pa]      : upstream pressure
% Pdo [pa]      : downstream pressure
% T_b  [K]      : upstream bulk temperature
% rho_b [kg/m3] : upstream bulk density
% c_v [J/kgK]   : specific heat at costant volume
% R_g [J/kgK]   : specific gas constant
% h_cl [m]      : clearance heigth
% L [m]         : clearance width
% 
% OUTPUT
% Mleak [kg/s] : leakage mass flow rate through path 3
%  
% HP: - Perfect gas
%     - Comnpressible flow
%     - Inviscid flow
%     - Isentropic nozzle
%
% Model inspired by: "Multi-Vane Expanders: Internal-Leakage Losses", O. Badr, 1985
% See Franzetti-Persico 2019 thesis for more info

%% DEFINITIONS %%
Cd    = 1;        % discharge coefficient
c_p   = c_v+R_g;  % gas specific heat at costant pressure [J/kgK]
Gamma = c_p/c_v;  % specific heat ratio [-]

% Correction in case of conditions in which Pdo is higher than Pup
% (overexpansion or overcompression). 
fact        = Pdo > Pup;
Pstar       = Pup(fact);
Pup(fact)   = Pdo(fact);
Pdo(fact)   = Pstar;

Pt = Pup*(2/(Gamma+1))^(Gamma/(Gamma-1));                 % critycal pressure in case of chocked flow
Pdo(Pdo./Pup < Pt./Pup) = Pt(Pdo./Pup < Pt./Pup);         % definition of the outlet pressure

Mleak = Cd.*h_cl.*W_cl.*rho_b.*sqrt(2*c_p*T_b.*((Pdo./Pup).^(2./Gamma)-(Pdo./Pup).^((Gamma+1)./Gamma))); % leakage mass flow rate through path 3 [kg/s]
Mleak(fact) = - Mleak(fact);
end

function [Mleak] = Ishii (Pup,Pdo,rho_b,mu_b,L_cl,h_cl,W_cl)
%This function compute leakages through axial clearance between vane side
%and end plates (path 3), using Ishii aproach.
% 
% INPUT
% Pup [pa]     : upstream pressure
% Pdo [pa]     : downstream pressure
% rho_b [kg/m^3]       : upstream density
% mu_b [Pa s]  : dynamic viscosity
% L_cl [m]     : clearance length
% h_cl [m]     : clearance height
% W_cl [m]     : clearance width
% 
% OUTPUT
% Mleak [kg/s] : leakage mass flow rate through path 3
%  
% HP: - Uncompresible flow
%     - Turbolent flow
%     - Viscous flow
%
% Model inspired by: "Refrigerant Leakage Flow Evaluation for Scroll Compressors”", N. Ishii, 1996
% See Franzetti-Persico 2019 thesis for more info
    
    %% CALCULATIONS %%
    % Initialization
    Mleak = NaN(length(Pup),1);
    x0 = 60;
    options.Display = 'off';
    
    % Correction in case of conditions in which Pdo is higher than Pup
    % (overexpansion or overcompression). 
    fact        = Pdo > Pup;
    Pstar       = Pup(fact);
    Pup(fact)   = Pdo(fact);
    Pdo(fact)   = Pstar;
    
    % Model
    for i = 1:length(Pup)
        Deq = h_cl*W_cl(i)./(2.*(h_cl + W_cl(i)));  % equivalent diameter
        FUN = @(Uo) 0.35*(4*Deq.*Uo.*rho_b(i)./mu_b(i)).^(-0.25).*L_cl.*Uo.^2./(8.*Deq) - (Pup(i)-Pdo(i))./rho_b(i);
        Uo  = fsolve(FUN, x0,options);
        Mleak(i) = Uo.*h_cl.*W_cl(i).*rho_b(i); % leakage mass flow rate through path 3 [kg/s]
    end
    Mleak(fact) = - Mleak(fact);
    
end

function [Mleak] = Suefuji(Pup,Pdo,T_b,h_cl,W_cl,mu_b,R_g,Gamma)  
%This function compute leakages through axial clearance between vane side
%and end plates (path 3), using Suefuji aproach.
%
% INPUT
% Pup [pa]     : upstream pressure
% Pdo [pa]     : downstream pressure
% T_b  [K]      : upstream bulk temperature
% h_cl [m]     : clearance height
% W_cl [m]     : clearance width
% mu_b [Pa s]  : dynamic viscosity
% R_g [J/kgK]   : specific gas constant
% Gamma [-]      : specific heat ratio
%
% OUTPUT
% Mleak [kg/s] : leakage mass flow rate through path 3
%
% HP: - Compresible flow
%     - Perfect gas
%     - Nozzle model with Friction
%     - Turbolent flow (laminar excluded for instability)
%        
% Model inspired by: "Performance Analysis of Hermetic Scroll Compressors", K. Suefuji  et Al., 1992 

    %% CALCULATIONS %%
    % Initialization
    Mleak = NaN(length(Pup),1);
    x0 = 0.001;
    options.Display = 'off';
    
    % Correction in case of conditions in which Pdo is higher than Pup
    % (overexpansion or overcompression). 
    fact        = Pdo > Pup;
    Pstar       = Pup(fact);
    Pup(fact)   = Pdo(fact);
    Pdo(fact)   = Pstar;
        
    for i=1:length(Pup)
        if (Pup(i)-Pdo(i))==0
            Mleak(i) = 0; % leakage mass flow rate through path 3 [kg/s]
        else
            Mleak(i) = fsolve( @(x)Suefuji_System(x,Pup(i),Pdo(i),T_b(i),h_cl,W_cl(i),mu_b(i),R_g,Gamma), x0 ,options); % leakage mass flow rate through path 3 [kg/s]
        end
    end
    Mleak(fact) = - Mleak(fact);
    
end

function FUN = Suefuji_System(x,Pup,Pdo,T_b,h_cl,W_cl,mu_b,R_g,Gamma)
    % Reynolds number [-]
    Re     = 2*x/(mu_b*(W_cl+h_cl));
    % Darcy factor [-] (Hy: only turbolent)
    lambda = 0.3164/Re^0.25;
    % polytropic index [-]
    n = Gamma*(1+lambda)/(1+Gamma*lambda);   
    % friction loss coefficient [-]
    csi = sqrt(((Pdo/Pup)^(2/n)-(Pdo/Pup)^((n+1)/n))/((Pdo/Pup)^(2/Gamma)-(Pdo/Pup)^((Gamma+1)/Gamma))); 
    
    FUN = x - csi*h_cl*W_cl*Pup/sqrt(R_g*T_b)*sqrt(2*Gamma/(Gamma-1)*((Pdo/Pup)^(2/Gamma)-(Pdo/Pup)^((Gamma+1)/Gamma)));

end

function [Mleak] = Yanagisawa (Pup,Pdo,T_b,c_v,L_cl,h_cl,W_cl,mu_b,R_g)
    %This function compute leakages through axial clearance between vane side
    %and end plates (path 3), using Yanagisawa & Shimizu approach
    % 
    % INPUT
    % Pup [pa]     : upstream pressure
    % Pdo [pa]     : downstream pressure
    % T_b  [K]      : upstream bulk temperature
    % c_v [J/kg K] : specific heat at costant volume
    % L_cl [m]       : clearance length
    % h_cl [m]     : clearance height
    % W_cl [m]     : clearance width
    % mu_b [Pa s]  : dynamic viscosity   
    % R_g [J/kgK]    : specific gas constant
    % 
    % OUTPUT
    % Mleak [kg/s] : leakage mass flow rate through path 3
    % 
    % NOTES: - subscripts stand for:
    %              up : upstream of clearance
    %              do : downstream of clearance
    %              t  : throat equivalernt channel
    %              e  : equivalent channel exit
    % HP: - Compresible flow
    %     - Ideal gas
    %     - Adiabatic Flow
    %     - Viscous flow, turbolent (laminar excluded for instability)
    %     - Fanno flow (entropy always increase. Max Mach is 1)
    %
    % Model inspired by: "Leakage losses with a rolling piston type compressor.
    % I: radial clearanceon the rolling pèiston", Yanagisawa & Shimizu, 1985
    %
    % See Franzetti-Persico 2019 for more info

    %% DEFINITIONS
    c_p   = c_v+R_g;  % gas specific heat at costant pressure [J/kgK]
    Gamma = c_p/c_v;  % specific heat ratio [-]
    options.Display = 'off';
    options.MaxFunctionEvaluations = 1e5;
    options.ConstraintTolerance = 1e-10;
    options.FunctionTolerance = 1e-10;
    options.StepTolerance     = 1e-10;
    options.Algorithm         = 'levenberg-marquardt';

    %% FIRST ATTEMPT SOLUTION: CHOKE HYPOTHESYS %%
    % Mach number at exit is set to 1. Then it is verified if the hypothesis is confirmed
    Mleak = NaN(length(Pup),1);
    x_guess = [0.2, 0.15]';   % initial guess
    
    for i=1:length(Pup)
        x = fsolve(@(x) Yanagawa_choke(x,T_b(i),W_cl(i),h_cl,R_g,mu_b(i),Gamma,L_cl,Pup(i)),x_guess,options); % solve system
        
        % solution estraction
        Mt    = x(1);
        Mleak(i) = x(2);
        
        % solve pressure ratios
        Pt_e   = (1/Mt)*((Gamma+1)/(2+(Gamma-1)*Mt.^2))^0.5;
        Pup_t  = (1+(Gamma-1)*Mt^2/2)^(Gamma/(Gamma-1));
        Pratio = Pup_t*Pt_e;
        
        %% COMPLETE SOLUTION %%
        % if condition below is verified, flow is NOT choking, and the solution
        % above is not valid. A new solution must be calculated.
        if Pratio > (Pup(i)/Pdo(i))
            x_guess = [Mt, 0.8, 0.001]'; % initial guess
            x  = fsolve(@(x) Yanagawa_NOchoke(x,T_b(i),W_cl(i),h_cl,R_g,mu_b(i),Gamma,L_cl,Pup(i),Pdo(i)),x_guess,options); % solve system
            
            % solution estraction
            Mleak(i) = x(3);
        end
    end
end
function FUN = Yanagawa_choke(x,T_b,W_cl,h_cl,R_g,mu_b,Gamma,L_cl,Pup)
    % This function set the system of 2 non linear equation under choking hypothesys
    % Unknowns of the system:
    % Mt = x(1)            : throat mach number [-] 
    % Mleak = x(2)     : leakage mass flow rate [kg/s]

    %% PREAMBLE %%
    FUN = NaN(2,1);
    Me  = 1;        

    % define lambda (Hy: turbolent flow)
    Re     = 2*x(2)/(mu_b*(W_cl+h_cl));
    lambda = 0.3164/Re^0.25;

    %% SYSTEM DEFINITION %%
    Deq     = 4*W_cl*h_cl/(2*(W_cl+h_cl));
    FUN (1) = (((1-x(1)^2)/(Gamma*x(1)^2)+(Gamma+1)/(2*Gamma)*log((Gamma+1)*x(1)^2/(2+(Gamma-1)*x(1)^2)))/(L_cl*lambda/Deq))-1; % I equation

    % Pressure ratios
    Pt_e   = (1/x(1))*((Gamma+1)/(2+(Gamma-1)*x(1).^2))^0.5; % ratio between throat pressure and exit pressure [-]
    Pup_t  = (1+(Gamma-1)*x(1)^2/2)^(Gamma/(Gamma-1));       % ratio between upstream pressure and throat pressure [-]
    Pratio = Pup_t*Pt_e;                                     % ratio between upstream and exit pressure [-]

    % Define exit condition
    Te = T_b/(1+Me^2*(Gamma-1)/2);  % temperature at channel exit [K]
    Pe = Pup/Pratio;                % pressure  at channel exit [Pa]
    Ve = Me*(Gamma*R_g*Te)^0.5;     % velocity at channel exit [m/s]

    % Mass flow rate 
    Mleak  = h_cl*W_cl*Ve*Pe/(R_g*Te);  % leakage mass flow rate [kg/s]
    FUN(2) = Mleak-x(2);                % II equation
end


function FUN = Yanagawa_NOchoke(x,T_b,W_cl,h_cl,R_g,mu_b,Gamma,L_cl,Pup,Pdo)
    % This function set the system of 3 non linear equation for not-choking flow
    % Unknowns of the system:
    % Mt = x(1)         : throat mach number [-]
    % Me = x(2)         : equivalent channel exit mach number [-]
    % Mleak = x(3)      : leakage mass flow rate [kg/s]

    %% PREAMBLE %%
    FUN = NaN(3,1);

     % define lambda (Hy: turbolent flow)
    Re     = 2*x(3)/(mu_b*(W_cl+h_cl));
    lambda = 0.3164/Re^0.25;

    %% SYSTEM DEFINITION %%
    Deq    = 4*W_cl*h_cl/(2*(W_cl+h_cl));
    Lf_i   = Deq/lambda*((1-x(1)^2)/(Gamma*x(1)^2)+(Gamma+1)/(2*Gamma)*log((Gamma+1)*x(1)^2/(2+(Gamma-1)*x(1)^2))); % new specific length [m]
    FUN(1) = (1-x(2)^2)/(Gamma*x(2)^2)+(Gamma+1)/(2*Gamma)*log((Gamma+1)*x(2)^2/(2+(Gamma-1)*x(2)^2))-(Lf_i-L_cl)*lambda/Deq;  % I equation

    % pressure ratios
    Pt_x   = (1/x(1))*((Gamma+1)/(2+(Gamma-1)*x(1).^2))^0.5; % ratio between throat pressure and local pressure [-]
    Pup_t  = (1+(Gamma-1)/2*x(1)^2)^(Gamma/(Gamma-1));       % ratio between upstream pressure and throat pressure [-]
    Pe_x   = (1/x(2))*((Gamma+1)/(2+(Gamma-1)*x(2).^2))^0.5; % ratio between exit pressure and local pressure [-]
    Pratio = Pup_t*Pt_x/Pe_x;                                % ratio between upstream and exit pressure [-]
    FUN(2) = abs(Pratio/(Pup/Pdo))-1;  % II equation

    % Define exit condition
    Te = T_b/(1+x(2)^2*(Gamma-1)/2); % temperature at channel exit [K]
    Pe = Pup/Pratio;                 % pressure  at channel exit [Pa]
    Ve = x(2)*(Gamma*R_g*Te)^0.5;    % velocity at channel exit [m/s]

    % Mass flow rate 
    Mleak  = h_cl*W_cl*Ve*Pe/(R_g*Te);  % leakage mass flow rate [kg/s]
    FUN(3) = abs(Mleak /x(3))-1;      % III equation
end