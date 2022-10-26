function [Fleak_g_PE_dis,Fleak_o_PE_dis,Fleak_g_PE_suc,Fleak_o_PE_suc,Fleak_g_PE_cell_closed,Fleak_o_PE_net,Fleak_g_PE_net,Dmg_PE_comp,FlowLeak_gas] = S3_LeakagePE(MUmodel_PE,PEmodel,process,d,d_hub,PEclr,dt_vano,pos_SucOpen,pos_SucClose,pos_DisOpen,pos_DisClose,Gamma,theta_vane,V_cell,V_g,m_g,m_o,rpm,p_suc,p_comp,p_del,T_g,~,R_g,mu_g_PE,c_v,Z,rho_o_PE,mu_o_PE,Npt_cell,Npt_comp,dt,s)
% This function compute leakages through clearance between rotor and end-plates (path 2)
%
% INPUT
% MUmodel_PE [string]  : viscosity model for two-phase gas-liquid mixture
% PEmodel [string]     : rotor-end-plate leakage model
% process              : process selected (1: compression  2: expansion
% d [m]                : rotor diameter
% d_hub [m]            : hub diameter
% PEclr [m]            : rotor-end-plate clearance size
% dt_vano [s]          : time needed for a vane to sweep a cell [s]
% pos_SucOpen [-]      : suction opening discretized angular position in vector theta
% pos_SucClose [-]     : suction closing discretized angular position in vector theta
% pos_DisOpen [-]      : delivery opening discretized angular position in vector theta
% pos_DisClose [-]     : delivery opening discretized angular position in vector theta
% Gamma [rad]          : cell angular aperture
% theta_vane [rad]     : array of discretized angular position  [0:2pi]
% V_cell [m^3]         : cell volumes array of a complete vane revolution [-gamma:2pi]
% V_g    [m^3]         : gas volume for each discretization point
% m_g    [kg]          : gas mass for each discretization point
% m_o    [kg]          : oil mass for each discretization point
% rpm [rpm]            : rotor angular velocity
% p_suc [Pa]           : suction pressure
% p_comp [Pa]          : vector of pressures from suction closure to discharge aperture
% p_del  [Pa]          : delivery pressure
% T_g [K]              : temperature array
% T_mix [K]            : oil + gas adiabatic mixing temperature for each discretization step
% R_g [J/kgK]          : specific gas constant
% mu_g_PE [Pa s]       : leaked gas dynamic viscosity
% c_v [J/kg K]         : specific heat at costant volume
% Z   [-]              : compressibility factor 
% rho_o_PE [kg/m3]     : leaked oil density
% mu_o_PE [Pa s]       : leaked oil kinematic viscosity
% Npt_cell [-]         : number of discretization points between 2 vanes
% Npt_comp [-]         : discretization points of closed cell process
% n_van [-]            : number of vanes
% dt [s]               : time step 
% boundary_condition   : type of boundary condition for Badr model

%
% OUTPUT
% Fleak_g_PE_dis [kg/s]  : leaking gas mass flow rate in discharge through path 2
% Fleak_o_PE_dis [kg/s]  : leaking oil mass flow rate in discharge through path 2
% Fleak_g_PE_suc [kg/s]  : leaking gas mass flow rate in suction through path 2
% Fleak_o_PE_suc [kg/s]  : leaking oil mass flow rate in suction through path 2
% Fleak_g_PE_net [kg/s]  : net leaking gas mass flow rate along compression process through path 2 
% Fleak_o_PE_net [kg/s]  : net leaking oil mass flow rate along compression process through path 2 
% Dmg_PE_comp [kg]       : cumulated gas mass variation due to leakage through path 2 along compression process
% Fleak_g_PE_net_sum [kg/s]: net leaking gas mass flow rate along compression process through path 2 

% HISTORY
% V08_Franzetti_Persico: creation of this section, development of the leakage model
% V09_Mellone: addition of Poiselle-Couette and Ishii model for PE leakage path and review of the viscosity models 
% V10.1_DeFranco_Genoni_Gianoncelli: review of suction and delivery time (further information available in the thesis) 

  %% DEFINITIONS %%
    % GEOMETRY
    omega = rpm*2*pi/60;  % angular velocity [rad/s]
    r     = d/2;          % rotor radius [m]
    r_hub = d_hub/2;      % shaft radius [m]

    % Angles  suction & delivery are expressed in function of
    % angles=[-Gamma,2*pi], need to subtract Npt_cell 
    theta_SucOpen = theta_vane(pos_SucOpen-Npt_cell);       % suction open [rad]   
    theta_SucClose   = theta_vane(pos_SucClose-Npt_cell);   % suction close [rad]
    theta_DisOpen = theta_vane(pos_DisOpen-Npt_cell);       % delivery open [rad]
    theta_DisClose    = theta_vane(pos_DisClose-Npt_cell);  % delivery close [rad]
    ARC_n_1=theta_DisClose;  % used to be consistent with name of variables given by Badr
    
    angle_dis = theta_vane(end)-theta_vane(pos_DisOpen-2*Npt_cell);  % Angle in discharge swept by reference blade [rad] 
    angle_suc = theta_vane(pos_SucClose-Npt_cell);                   % Angle in suction swept by reference blade [rad]
    angle_cell_closed=theta_vane(pos_DisOpen-2*Npt_cell)-theta_vane(pos_SucClose-Npt_cell); % Angle in compression/expansion swept by reference blade [rad]
        
    dt_angle_dis  = angle_dis/omega;                        % A closed cell see discharge port for this amount of time [s]
    dt_angle_suc  = (angle_suc)/omega + dt_vano;            % A closed cell see suction port for this amount of time(include both phase I and phase II in suction, see below) [s]
    dt_angle_net=angle_cell_closed/omega;                   % A closed cell is completely sealed (closed cell phase) for this amount of time [s] 
    
    % THERMODYNAMICS
    
    % Pressure
    p_in = p_suc;            % suction pressure [Pa]

    if isnan(p_del)          % definition of delivery pressure
        p_out = p_comp(end);
    else
        p_out = p_del;
    end

    % Mass fraction
    X_g = m_g./(m_g+m_o);                                           % gas mass fraction in the current cell [-]

    % Volume fraction
    X_vol   = V_g./V_cell(pos_SucClose:pos_DisOpen-Npt_cell)';      % gas volume fraction in the current cell [-]
    X_vol_m = mean(X_vol);                                          % mean gas volume fraction [-]

    % Density
    rho_g   = p_comp./(Z.*R_g.*T_g);                                % gas density in the current cell [kg/m^3]
    rho_b   = 1./sum([X_g./rho_g,(1-X_g)./rho_o_PE],2,'omitnan');   % bulk density in the current cell [kg/m^3]
    rho_b_m = mean(rho_b);                                          % average density [kg/m3]
     
    % Viscosity model
    switch MUmodel_PE
        case 'Dukler'
            mu_b  = rho_b.*sum([X_g.*mu_g_PE./rho_b,(1-X_g).*mu_o_PE./rho_o_PE],2,'omitnan');            
        case 'Maxwell'
            mu_b_M1 = mu_o_PE.*sum([2*mu_o_PE*ones(Npt_comp,1),mu_g_PE-2.*X_g.*(mu_o_PE-mu_g_PE)],2,'omitnan')./sum([2*mu_o_PE*ones(Npt_comp,1),mu_g_PE+X_g.*(mu_o_PE-mu_g_PE)],2,'omitnan');          % -Maxwell 1 (gas based)
            mu_b_M2 = mu_g_PE.*sum([2*mu_g_PE,mu_o_PE-2.*(1-X_g).*(mu_g_PE-mu_o_PE)],2,'omitnan')./sum([2*mu_g_PE,mu_o_PE+(1-X_g).*(mu_g_PE-mu_o_PE)],2,'omitnan');  % -Maxwel 2 (liquid based)
            mu_b    = sum([mu_b_M1,mu_b_M2],2,'omitnan')/2;  
            if sum(isnan(mu_b_M1)) >= 1
                mu_b =  mu_b_M2;
            end
        case 'Awad'
            mu_b    = 1/4*(sum([(3*X_g-1).*mu_g_PE,(3*(1-X_g)-1).*mu_o_PE],2,'omitnan') + sqrt(sum([sum([(3*X_g-1).*mu_g_PE,(3*(1-X_g)-1).*mu_o_PE],2,'omitnan').^2,8*mu_o_PE.*mu_g_PE],2,'omitnan')));
    end
    
    mu_b_m = mean(mu_b);  % dynamic viscosity [Pa s]   

   switch PEmodel
       
       case 'Badr'
           [FlowLeak,FlowLeak_zero] = Badr_model(process,d,d_hub,PEclr,Gamma,theta_vane,V_cell,p_comp,p_del,R_g,c_v,Npt_cell,mu_b_m, rho_b_m,p_out,p_in,omega,r,r_hub,theta_SucOpen,theta_SucClose,theta_DisOpen,theta_DisClose,ARC_n_1);
       case 'Poiselle-Couette'
            [FlowLeak,FlowLeak_zero] = Poiselle_Couette(process,d,d_hub,PEclr,Gamma,theta_vane,V_cell,p_comp,p_del,R_g,c_v,Npt_cell,mu_b_m, rho_b_m,p_out,p_in,omega,r,r_hub,theta_SucOpen,theta_SucClose,theta_DisOpen,theta_DisClose,ARC_n_1,pos_SucClose,pos_SucOpen,pos_DisClose,pos_DisOpen,mu_b,rho_b,s,X_g,rho_o_PE,rho_g);
       case 'Ishii'
            [FlowLeak,FlowLeak_zero] = Ishii_model(process,d,d_hub,PEclr,Gamma,theta_vane,V_cell,p_comp,p_del,R_g,c_v,Npt_cell,mu_b_m, rho_b_m,p_out,p_in,omega,r,r_hub,theta_SucOpen,theta_SucClose,theta_DisOpen,theta_DisClose,ARC_n_1,pos_SucClose,pos_SucOpen,pos_DisClose,pos_DisOpen,mu_b,rho_b,s,X_g);
            
   end
   
     %% COMPUTATION %%
     % Leakage contribution from gas and oil through path 2 for each discretization in [0;2*pi] (considering a mean volume fraction)  [kg/s] 
       FlowLeak_gas = FlowLeak.*X_vol_m;
       FlowLeak_oil = FlowLeak.*(1-X_vol_m); 
       
       %leakage contribution considering volume reduction in last cell, to
       %be used to calculate effective contribution in suction and
       %discharge of the leakages through path 2
       FlowLeak_zero_gas=FlowLeak_zero.*X_vol_m;
       FlowLeak_zero_oil = FlowLeak_zero.*(1-X_vol_m);  
     
       Fleak_g_PE_net = (FlowLeak_gas(pos_SucClose-Npt_cell:pos_DisOpen-2*Npt_cell))';  % net leaking gas mass flow rate along compression process through path 2 [kg/s] 
       Fleak_o_PE_net = (FlowLeak_oil(pos_SucClose-Npt_cell:pos_DisOpen-2*Npt_cell))';  % net leaking oil mass flow rate along compression process through path 2 [kg/s]

       Dmg_PE_comp    = cumsum(Fleak_g_PE_net*dt); % cumulated gas mass variation due to leakage through path 2 along compression process [kg]
     % Dmo_VS_comp  = cumsum(Fleak_o_PE_net*dt);     cumulated oil mass variation due to leakage through path 2 along compression process [kg]
               
          switch process
              
              case 1   % compression process
                  

          Fleak_g_PE_suc_I=FlowLeak_gas-FlowLeak_zero_gas;  %difference between leakages considering whole cell volume in last part of compression/discharge and leakages considering volume reduction of cell
          Fleak_g_PE_suc_I=Fleak_g_PE_suc_I(end-Npt_cell:end);
          Fleak_g_PE_suc_II = FlowLeak_gas(1:pos_SucClose-Npt_cell-1);
          Fleak_g_PE_suc    = (sum([Fleak_g_PE_suc_I,Fleak_g_PE_suc_II]*dt)/dt_angle_suc)';   % leaking gas mass flow rate in suction through path 2 [kg/s] 
           
          Fleak_o_PE_suc_I=FlowLeak_oil-FlowLeak_zero_oil;  %difference between leakages considering whole cell volume in last part of compression/discharge and leakages considering volume reduction of cell
          Fleak_o_PE_suc_I=Fleak_o_PE_suc_I(Fleak_o_PE_suc_I>0);
          Fleak_o_PE_suc_II = FlowLeak_oil(1:pos_SucClose-Npt_cell-1);
          Fleak_o_PE_suc    = (sum([Fleak_o_PE_suc_I,Fleak_o_PE_suc_II]*dt)/(dt_angle_suc))';   % leaking oil mass flow rate in suction through path 2 [kg/s]
             
          Fleak_g_PE_dis = (sum([FlowLeak_gas(pos_DisOpen-2*Npt_cell+1:end-Npt_cell),FlowLeak_zero_gas(end-Npt_cell:end)]*dt)/dt_angle_dis)'; % leaking gas mass flow rate in discharge [kg/s]
          Fleak_o_PE_dis = (sum([FlowLeak_oil(pos_DisOpen-2*Npt_cell+1:end-Npt_cell),-Fleak_o_PE_suc_I(end:-1:1)]*dt)/dt_angle_dis)'; % leaking gas mass flow rate in discharge [kg/s]

      
              case 2      % expansion process
              
          Fleak_g_PE_suc_I=FlowLeak_gas-FlowLeak_zero_gas; %difference between leakages considering whole cell volume in last part of compression/discharge and leakages considering volume reduction of cell
          Fleak_g_PE_suc_I=Fleak_g_PE_suc_I(end-Npt_cell:end);
          Fleak_g_PE_suc_II = FlowLeak_gas(1:pos_SucClose-Npt_cell-1);
          Fleak_g_PE_suc    = (sum([Fleak_g_PE_suc_I,Fleak_g_PE_suc_II]*dt)/dt_angle_suc)';   % leaking gas mass flow rate in suction through path 2 [kg/s] 
           
          Fleak_o_PE_suc_I=FlowLeak_oil-FlowLeak_zero_oil;  %difference between leakages considering whole cell volume in last part of compression/discharge and leakages considering volume reduction of cell
          Fleak_o_PE_suc_I=Fleak_o_PE_suc_I(Fleak_o_PE_suc_I>0);
          Fleak_o_PE_suc_II = FlowLeak_oil(1:pos_SucClose-Npt_cell-1);
          Fleak_o_PE_suc    = (sum([Fleak_o_PE_suc_I,Fleak_o_PE_suc_II]*dt)/(dt_angle_suc))';   % leaking oil mass flow rate in suction through path 2 [kg/s]
             
          Fleak_g_PE_dis = (sum([FlowLeak_gas(pos_DisOpen-2*Npt_cell+1:end-Npt_cell),FlowLeak_zero_gas(end-Npt_cell:end)]*dt)/dt_angle_dis)'; % leaking gas mass flow rate in discharge [kg/s]
          Fleak_o_PE_dis = (sum([FlowLeak_oil(pos_DisOpen-2*Npt_cell+1:end-Npt_cell),-Fleak_o_PE_suc_I(end:-1:1)]*dt)/dt_angle_dis)'; % leaking gas mass flow rate in discharge [kg/s]

      
          end
          
 
    Fleak_g_PE_cell_closed=Dmg_PE_comp(end)/(dt_angle_net); % net leaking gas mass flow rate along compression process through path 2 [kg/s] 
    
    % removing the last element since correspond to the opening of the delivery
    Fleak_g_PE_net(end) = [];
    Fleak_o_PE_net(end) = [];
    
    
function [FlowLeak,FlowLeak_zero] = Badr_model(process,~,~,PEclr,Gamma,theta_vane,V_cell,p_comp,~,R_g,c_v,Npt_cell,mu_b_m, rho_b_m,p_out,p_in,omega,r,r_hub,theta_SucOpen,theta_SucClose,theta_DisOpen,theta_DisClose,ARC_n_1)
%This function compute leakages through axial clearances between the rotor and the two end- plates using Badr compressible approach. Model inspired by:
%"Multi-Vane Expanders: Internal-Leakage Losses", O. Badr, 1985 (see Franzetti-Persico thesis and Mellone thesis for more info)
%Three different boundary condition are proposed: linear, adiabatic, effective 


          %boundary condition setting
          boundary_condition='effective';          %new boundary condition for BADR formula: linear, adiabatic, effective
       
          
          % Dimensionless number
          Re_PEclr = rho_b_m*omega*PEclr^2/mu_b_m;           % Reynolds number
          Eu       = (p_in-p_out)/(omega^2*rho_b_m*r);       % Euler number
          % S        = r/PEclr;                              % Gap aspect ratio
          beta     = r_hub/r;                                % shaft-rotor radius ratio
          
   %% SOLVING EQUATION
          % fourier coefficients calculate with coeff_Fourier function base on use of 
          % symbolic tools, four boundary condition avaiable: linear,
          % adiabatic, effective

          % number of harmonically related sinusoids of Fourier series
          N=50;

    % BOUNDARY CONDITION at r* = 1 (rotor outlet radius)
    
    switch boundary_condition
 
          case 'linear'
              
          % LINEAR compressor resampled with reference to trailing vane,
          % zero at tangency position: BADR
         
          % Pressure trend at rotor surface (hypotized to be linear)
          % A0, AK, Bk have been calculated by fourier seriers using the
          % following expression of boundary condition
          %{
          
           % Initialization
          f_linear = NaN(size(theta_vane));
          
          fact1 = theta_vane <= theta_SucClose;
          fact2 = theta_vane > theta_SucClose & theta_vane <= (theta_DisOpen-Gamma);
          fact3 = theta_vane > (theta_DisOpen-Gamma) & theta_vane <= (theta_DisClose-Gamma);
          fact4 = theta_vane > (theta_DisClose-Gamma) & theta_vane <= ARC_n_1;
          fact5= theta_vane> ARC_n_1;
          
          f_linear(fact1) = 0;
          f_linear(fact2) = (-theta_SucClose+theta_vane(fact2))./(theta_DisOpen-theta_SucClose-Gamma);
          f_linear(fact3) = 1;
          f_linear(fact4) = (-theta_vane(fact4)+ARC_n_1)./(ARC_n_1-theta_DisClose+Gamma);
          f_linear(fact5)=0;
          %}
          
          k=1:N;
         
          Ak=1./(pi.*(1+beta.^(2.*k))).*((sin(k.*(Gamma - theta_DisOpen)) - sin(k.*(Gamma - theta_DisClose)))./k + (sin(k.*(theta_SucOpen - Gamma + 2.*pi)) + sin(k.*(Gamma - theta_DisClose)))./k - (sin(k.*(theta_SucOpen - Gamma + 2.*pi)) - sin(2.*pi.*k))./k - (cos(k.*(Gamma - theta_DisOpen)) - cos(k.*theta_SucClose) + k.*(Gamma.*sin(k.*(Gamma - theta_DisOpen)) - theta_DisOpen.*sin(k.*(Gamma - theta_DisOpen)) + theta_SucClose*sin(k.*(Gamma - theta_DisOpen))))./(k.^2.*(Gamma - theta_DisOpen + theta_SucClose)));
          Bk=1./(pi.*(1+beta.^(2.*k))).*((cos(k.*(Gamma - theta_DisOpen)) - cos(k.*(Gamma - theta_DisClose)))./k + (cos(k.*(theta_SucOpen - Gamma + 2.*pi)) - cos(2.*pi.*k))./k - (cos(k.*(theta_SucOpen - Gamma + 2.*pi)) - cos(k.*(Gamma - theta_DisClose)))./k + (sin(k.*theta_SucClose) + sin(k.*(Gamma - theta_DisOpen)) - k.*(Gamma*cos(k.*(Gamma - theta_DisOpen)) - theta_DisOpen.*cos(k.*(Gamma - theta_DisOpen)) + theta_SucClose.*cos(k.*(Gamma - theta_DisOpen))))./(k.^2.*(Gamma - theta_DisOpen + theta_SucClose)));
 
          FourArg = (1-beta.^(2*k')).*(Ak'.*(sin(k'*theta_vane+Gamma/Npt_cell)-sin(k'*(theta_vane)))-Bk'.*(cos(k'*theta_vane+Gamma/Npt_cell)-cos(k'*(theta_vane))));
       
         % Leakage mass flow rate through path 2 [kg/s] for each
         % discretization step (trailing vane as reference and tangency as
         % starting point)
         
          FlowLeak_punctual =  process*Eu*Re_PEclr*omega*r^2*PEclr*rho_b_m/6*sum(FourArg);
         
        
          FlowLeak_complete=[FlowLeak_punctual,FlowLeak_punctual(1:Npt_cell)];  %leakages along whole compression proces from 0 to 2pi+gamma
          FlowLeak = movsum(FlowLeak_complete,[0 Npt_cell],'Endpoints','discard');
          
          % Corerection to account for the numerical error linked to the usage of Fourier series instead of the analytical solution
         
          eps_corr = sum(FlowLeak)/length(FlowLeak);
          FlowLeak = FlowLeak - eps_corr;
          
          %calculation of leakage must be repeated in order to calculate
          %leakage in last cell considering the reduction of volume, so
          %only the contribute of discharge part and not the part of the
          %cell that will be in suction
          FlowLeak_zero= movsum(FlowLeak_punctual,[0 Npt_cell],'Endpoints','shrink');
       
     
          case 'adiabatic'
              
          % Adiabatic resampled with reference to trailing vane,
          % zero at tangency position
    
          % Initialization
          
          f_adiabatic = NaN(size(theta_vane));
          k_adiabatic               = 1+R_g/c_v; 
          
          % Pressure trend at rotor surface (hypotized to be linear) 
          
          fact1 = theta_vane < theta_SucClose;
          fact2 = theta_vane >= theta_SucClose & theta_vane <= (theta_DisOpen-Gamma);
          fact3 = theta_vane > (theta_DisOpen-Gamma) & theta_vane <= (theta_DisClose-Gamma);
          fact4 = theta_vane > (theta_DisClose-Gamma) & theta_vane <= ARC_n_1;
          fact5= theta_vane> ARC_n_1;
    
          V_ad=V_cell(Npt_cell:end);  %volume of compression/expansion chamber along whole process
          Volume=V_ad(fact2);         %volume of compression/expansion chamber during closed cell phase
          
          switch process
              
              case 1
                  
          f_adiabatic(fact1) = 0;
          f_adiabatic(fact2) = (((Volume(1)./Volume).^k_adiabatic)-p_in./10.^5)./(p_out./10^5);
          f_adiabatic(fact3) = 1;
          f_adiabatic(fact4) = 1;
          f_adiabatic(fact5) = 1;
          
              case 2
                  
          f_adiabatic(fact1) = 1;
          f_adiabatic(fact2) = (((Volume(end)./Volume).^k_adiabatic)-p_out./10.^5)./(p_in./10^5);
          f_adiabatic(fact3) = 0;
          f_adiabatic(fact4) = (-(theta_DisClose-Gamma)+theta_vane(fact4))./(ARC_n_1-(theta_DisClose-Gamma));
          f_adiabatic(fact5) = 1; 
          
          end
          
          k=1:N;
                
          Ak=(trapz(theta_vane,f_adiabatic.*(cos(k'*theta_vane)),2))'./(pi.*(1+beta.^(2*k)));
          Bk=(trapz(theta_vane,f_adiabatic.*(sin(k'.*theta_vane)),2))'./(pi.*(1+beta.^(2.*k)));
         
         FourArg = (1-beta.^(2*k')).*(Ak'.*(sin(k'*(theta_vane+Gamma/Npt_cell))-sin(k'*(theta_vane)))-Bk'.*(cos(k'*(theta_vane+Gamma/Npt_cell))-cos(k'*(theta_vane))));
        
         % Leakage mass flow rate through path 2 [kg/s] for each
         % discretization step (trailing vane as reference and tangency as
         % starting point)
         
          FlowLeak_punctual =  process*Eu*Re_PEclr*omega*r^2*PEclr*rho_b_m/6*sum(FourArg);
         
        
          FlowLeak_complete=[FlowLeak_punctual, FlowLeak_punctual(1:Npt_cell)];  %leakages along whole compression proces from 0 to 2pi+gamma
          FlowLeak = movsum(FlowLeak_complete,[0 Npt_cell],'Endpoints','discard');
          
          % Corerection to account for the numerical error linked to the usage of Fourier series instead of the analytical solution
         
          eps_corr = sum(FlowLeak)/length(FlowLeak);
          FlowLeak = FlowLeak - eps_corr;
          
          %calculation of leakage must be repeated in order to calculate
          %leakage in last cell considering the reduction of volume, so
          %only the contribute of discharge part and not the part of the
          %cell that will be in suction
          FlowLeak_zero= movsum(FlowLeak_punctual,[0 Npt_cell],'Endpoints','shrink');
         
     
          case 'effective'
              
          % Effective: fact 4 and fact 5 adapted to expected pressure trend
         
          %Initialization
          
          f_effective = NaN(size(theta_vane));
    
          fact1 = theta_vane < theta_SucClose;
          fact2 = theta_vane >= theta_SucClose & theta_vane <= (theta_DisOpen-Gamma);
          fact3 = theta_vane > (theta_DisOpen-Gamma) & theta_vane <= (theta_DisClose-Gamma);
          fact4 = theta_vane > (theta_DisClose-Gamma) & theta_vane <= ARC_n_1;
          fact5= theta_vane> ARC_n_1;
          
          switch process
              case 1

          f_effective(fact1)=0;
          f_effective(fact2)=(((p_comp(1:end))./10^5)-p_in./10.^5)./(p_out./10.^5);   %%%cambiare la PRESSIONE!!!!!!!!
          f_effective(fact3)=1;
          f_effective(fact4) = 1;
          f_effective(fact5) = 1;
          
              case 2

          f_effective(fact1)=1;
          f_effective(fact2)=(((p_comp(1:end))./10^5)-p_out./10.^5)./(p_in./10.^5);
          f_effective(fact3)=0;
          f_effective(fact4) =(-(theta_DisClose-Gamma)+theta_vane(fact4))./(ARC_n_1-(theta_DisClose-Gamma));
          f_effective(fact5) = 1;
          
          end
          
          N=50;
          k=1:N;

          Ak=(trapz(theta_vane,f_effective.*(cos(k'*theta_vane)),2))'./(pi.*(1+beta.^(2*k)));
          Bk=(trapz(theta_vane,f_effective.*(sin(k'.*theta_vane)),2))'./(pi.*(1+beta.^(2.*k)));
          
          FourArg = (1-beta.^(2*k')).*(Ak'.*(sin(k'*(theta_vane+Gamma/Npt_cell))-sin(k'*(theta_vane)))-Bk'.*(cos(k'*(theta_vane+Gamma/Npt_cell))-cos(k'*(theta_vane))));
         
         % Leakage mass flow rate through path 2 [kg/s] for each
         % discretization step (trailing vane as reference and tangency as
         % starting point)

          FlowLeak_punctual =  process*Eu*Re_PEclr*omega*r^2*PEclr*rho_b_m/6*sum(FourArg);
         
          FlowLeak_complete=[FlowLeak_punctual, FlowLeak_punctual(1:Npt_cell)];   %leakages along whole compression proces from 0 to 2pi+gamma
          FlowLeak = movsum(FlowLeak_complete,[0 Npt_cell],'Endpoints','discard');
          
          % Corerection to account for the numerical error linked to the usage of Fourier series instead of the analytical solution
          eps_corr = sum(FlowLeak)/length(FlowLeak);
          FlowLeak = FlowLeak - eps_corr;
          
          %calculation of leakage must be repeated in order to calculate
          %leakage in last cell considering the reduction of volume, so
          %only the contribute of discharge part and not the part of the
          %cell that will be in suction
          FlowLeak_zero= movsum(FlowLeak_punctual,[0 Npt_cell],'Endpoints','shrink');
      
    end
   
    
    %% CHECKS %%
    % Conformity of the model applied (laminar flow to solve Navier - Stokes equations
%     chck = Re_PEclr*S^2 > 1e5;
%     
%     if chck
%         warning( 'S3_LeakagePE:Badr','Leaking flow should not be considered laminar according to Daily and Nage. Badr is not the best model.');
%         SX_Logfile ('w',{lastwarn});
%     end
    
   %integral on the circumference must be equal to zero, conservation of
   %mass, the mass that enter in the cells must be equal to mass that exit
   %from cells
   
   mass_variation=sum(FlowLeak);
   mass_chck=abs(mass_variation)>1e-9;
   
   if mass_chck
   warning( 'S3_LeakagePE:Badr','Leaking flow does not respect the mass balance');
         SX_Logfile ('w',{lastwarn});
   end
   
end

function [FlowLeak,FlowLeak_zero] = Poiselle_Couette(process,~,~,PEclr,Gamma,theta_vane,~,p_comp,~,~,~,Npt_cell,~, ~,p_out,p_in,~,r,~,~,theta_SucClose,theta_DisOpen,theta_DisClose,ARC_n_1,pos_SucClose,~,~,pos_DisOpen,mu_b,rho_b,s,X_g,rho_o_PE,rho_g)
%This function compute leakages through axial clearances between the rotor and the two end- plates using Poiselle_Couette approach. Model inspired by:
%"Experimental and Theoretical Investigation of a Sliding Vane Compressor-Expander Unit for an R-134a Automotive Vapour Compression Refrigeration System", Chukwudi Azih, 2007
%(see Mellone thesis for more info)
%the circumference along which the leakages takes place is supposed to be
%divided in N streamtubes without any interaction between each
%other, it is assumed that all the leaking mass will go towards
%suction in compressor and discharge in expander

  switch process
        
        case 1   %compression process
            
         N_points=length(theta_vane)-(pos_SucClose-Npt_cell);                 % number of streamtubes, taken as number of points between suction close and end of discretization range
         N_points_short=pos_SucClose-Npt_cell;                                % number of discretization points in suction, it is different from N_points       
         distance_out=(theta_vane(end)-theta_SucClose)/(N_points-1);          % angle between each streamtubes in discharge and compression phase
         distance_in=(theta_SucClose)/(N_points-1);                           % angle between each streamtubes in suction phase
         theta_int=0:((theta_SucClose)/(N_points_short-1)):(theta_SucClose);  % discretization of suction phase according to N_points_short
         theta_2=(theta_SucClose):distance_out:theta_vane(end);               % discretization of discharge and compression phase according to number of streamtubes
         theta_1=0:distance_in:(theta_SucClose);                              % discretization of suction phase according to number of streamtubes
         L=2*r*sin((distance_out+distance_in)./4);                            % width of a streamtubes, takes as the mean between inlet and outlet
         
        case 2   %expansion process
            
         N_points=length(theta_vane)-(pos_DisOpen-2*Npt_cell);                          % number of streamtubes, taken as number of points between suction close and end of discretization range
         N_points_short=pos_DisOpen-2*Npt_cell;                                         % number of discretization points in suction, it is different from N_points
         distance_out=(theta_vane(end)-(theta_DisOpen-Gamma))/(N_points-1);             % angular distance between each streamtubes in discharge phase
         distance_in=(theta_DisOpen-Gamma)/(N_points-1);                                % angular distance between each streamtubes in suction and expansion phase
         theta_int=0:((theta_DisOpen-Gamma)/(N_points_short-1)):(theta_DisOpen-Gamma);  % discretization of discharge phase according to number of streamtubes
         theta_2=(theta_DisOpen-Gamma):distance_out:theta_vane(end);                    % discretization of discharge phase according to N_points_short
         theta_1=0:distance_in:(theta_DisOpen-Gamma);                                   % discretization of suction phase according to number of streamtubes
         L=2*r*sin((distance_out+distance_in)./4);                                      % width of a streamtubes, take as the mean between inlet and outlet
  end
    
   
    delta_theta=theta_2(end:-1:1)-theta_1; % [rad]
    fact_delta=delta_theta>180;      
    delta_theta(fact_delta)=2*pi-delta_theta(fact_delta); 
    W_PC=2.*r.*sin(delta_theta./2); %length of each streamtubes
    
    % check if distance is lower than vane thickness (s), that is considered
    % as the minimum length of a streamtube
    
    flag_W=W_PC<s;
    W_PC(flag_W)=s;
    
    %boundary condition, same approach of Badr effective model

          f_effective = NaN(size(theta_vane)); % initialization
    
          fact1 = theta_vane <= theta_SucClose;
          fact2 = theta_vane > theta_SucClose & theta_vane <= (theta_DisOpen-Gamma);
          fact3 = theta_vane > (theta_DisOpen-Gamma) & theta_vane <= (theta_vane(end)-Gamma);
          fact4 = theta_vane > (theta_vane(end)-Gamma) & theta_vane <= theta_vane(end);
          fact5= theta_vane> ARC_n_1;
          
          switch process
              
              case 1  % compression process
                  
          f_effective(fact1)=0;
          f_effective(fact2)=(((p_comp(length(p_comp)-sum(fact2)+1:end)))-p_in)./(p_out);   
          f_effective(fact3)=1;
          f_effective(fact5) = 1;
          f_effective(fact4)=1;
          
          P_PC=(f_effective((sum(fact1))+1:end))*(p_out-p_in); % trend of relative pressure in compression and discharge phases
          
          rho_b_pc=mean(rho_b);
          mu_b_pc=mean(mu_b);
          FlowLeak_PC=rho_b_pc.*L.*PEclr.^3.*(P_PC)./(12.*mu_b_pc.*W_PC(end:-1:1)); % leaking mass flow rate in compression and discharge phases
          
          FlowLeak_PC_resampled=interp1(theta_1,FlowLeak_PC,theta_int,'linear');  %FlowLeak_PC is calculated according to N_points, need to resample it according to length of theta_int (N_points_short) (punctual)
          FlowLeak_PC_resampled=FlowLeak_PC_resampled.*(sum(FlowLeak_PC)/sum(FlowLeak_PC_resampled));  %coorection of the results to respect the conservation of mass (punctual)
          
          FlowLeak_1=[-FlowLeak_PC_resampled(end:-1:1)];       % leaking mass flow rate along suction phases (punctual)
          FlowLeak_0=-FlowLeak_PC_resampled(Npt_cell:-1:1);    % leaking mass flow rate along last cell of suction phases  (punctual)
          FlowLeak_2=[FlowLeak_PC, FlowLeak_1(1:Npt_cell)]; % leaking mass flow rate along compression and discharge phases, considering also that part of the last cell will be divided between suction and discharge (punctual)
          
          FlowLeak_suction = -movsum(FlowLeak_1,[0 Npt_cell],'Endpoints','discard');   % leaking mass flow rate along suction phases
          FlowLeak_suction_end=-movsum(FlowLeak_0,[0 Npt_cell],'Endpoints','shrink');      % leaking mass flow rate along last cell of suction phases 
          FlowLeak_comp_dis = -movsum(FlowLeak_2,[0 Npt_cell],'Endpoints','discard');   % leaking mass flow rate along compression and discharge phases, considering also that part of the last cell will be divided between suction and discharge
          FlowLeak=[FlowLeak_suction,FlowLeak_suction_end,FlowLeak_comp_dis];                                % leaking mass flow rate along whole circumference
          
          % Corerection to account for the numerical error linked to the usage of Fourier series instead of the analytical solution
          eps_corr = sum(FlowLeak)/length(FlowLeak);
          FlowLeak = FlowLeak - eps_corr;
          
          %calculation of leakage must be repeated in order to calculate
          %leakage in last cell considering the reduction of volume, so
          %only the contribute of discharge part and not the part of the
          %cell that will be in suction
          FlowLeak_zero = [ FlowLeak_suction,FlowLeak_suction_end, -movsum(FlowLeak_PC,[0 Npt_cell],'Endpoints','shrink')];

              case 2   % expansion process
                  
          f_effective(fact1)=1;
          f_effective(fact2)=(((p_comp(length(p_comp)-sum(fact2)+1:end)))-p_out)./(p_in);   
          f_effective(fact3)=0;
          f_effective(fact4)=0;
          f_effective(fact5) = 0;
          
          p_high=(f_effective(1:(sum(fact2)+sum(fact1))))*((p_in)-p_out);       % trend of relative pressure in suction and expansion phase
          P_PC=interp1(theta_int,p_high,theta_1,'linear');                 % pressure resampled according to number of discretization points in suction
        
          rho_b_pc=mean(rho_b);
          mu_b_pc=mean(mu_b);
          
          FlowLeak_PC=rho_b_pc.*L.*PEclr.^3.*(P_PC)./(12.*mu_b_pc.*W_PC);  
          FlowLeak_PC_resampled=interp1(theta_1,FlowLeak_PC,theta_int,'linear');
          FlowLeak_PC_resampled=FlowLeak_PC_resampled.*(sum(FlowLeak_PC)/sum(FlowLeak_PC_resampled));
          
          FlowLeak_1=FlowLeak_PC;       % leaking mass flow rate along suction and expansion phases (punctual)
          FlowLeak_0=FlowLeak_PC(end-Npt_cell:1:end);    % leaking mass flow rate along last cell of expansion phases  (punctual)
          FlowLeak_2=[-FlowLeak_PC_resampled(end:-1:1), FlowLeak_1(1:Npt_cell)]; % leaking mass flow rate along c discharge phases, considering also that part of the last cell will be divided between suction and discharge (punctual)
          
          FlowLeak_suc_exp = -movsum(FlowLeak_1,[0 Npt_cell],'Endpoints','discard');   % leaking mass flow rate along suction and expansion phases
          FlowLeak_exp_end=-movsum(FlowLeak_0,[0 Npt_cell],'Endpoints','shrink');      % leaking mass flow rate along last cell of expansion phases  
          FlowLeak_discharge = -movsum(FlowLeak_2,[0 Npt_cell],'Endpoints','discard');   % leaking mass flow rate along c discharge phases, considering also that part of the last cell will be divided between suction and discharge
          FlowLeak=[FlowLeak_suc_exp,FlowLeak_exp_end(2:end),FlowLeak_discharge];                                % leaking mass flow rate along whole circumference
          
          % Corerection to account for the numerical error linked to the usage of Fourier series instead of the analytical solution
          eps_corr = sum(FlowLeak)/length(FlowLeak);
          FlowLeak = FlowLeak - eps_corr;
          
          %calculation of leakage must be repeated in order to calculate
          %leakage in last cell considering the reduction of volume, so
          %only the contribute of discharge part and not the part of the
          %cell that will be in suction
          FlowLeak_zero = [ FlowLeak_suc_exp,FlowLeak_exp_end(2:end), movsum(FlowLeak_PC_resampled(end:-1:1),[0 Npt_cell],'Endpoints','shrink')];

          end
    %integral on the circumference must be equal to zero, conservation of
   %mass, the mass that enter in the cells must be equal to mass that exit
   %from cells
   
   mass_variation=sum(FlowLeak);
   mass_chck=abs(mass_variation)>1e-9;
   
   if mass_chck
   warning( 'S3_LeakagePE:Poiselle_Couette','Leaking flow does not respect the mass balance');
         SX_Logfile ('w',{lastwarn});
   end
 
   
end

function [FlowLeak,FlowLeak_zero] = Ishii_model(process,~,~,PEclr,Gamma,theta_vane,~,p_comp,~,~,~,Npt_cell,~, ~,p_out,p_in,~,r,~,~,theta_SucClose,theta_DisOpen,theta_DisClose,ARC_n_1,pos_SucClose,~,~,pos_DisOpen,mu_b,rho_b,s,X_g)
%This function compute leakages through axial clearances between the rotor and the two end- plates using Ishii approach. Model inspired by:
%"Refrigerant Leakage Flow Evaluation for Scroll Compressors", N. Ishii, 1996 (see Mellone thesis for more info)
%the circumference along which the leakages takes place is supposed to be divided in N streamtubes without any interaction between each
%other, it is assumed that all the leaking mass will go towards suction   

    switch process
        
          case 1   %compression process
            
         N_points=length(theta_vane)-(pos_SucClose-Npt_cell);                 % number of streamtubes, taken as number of points between suction close and end of discretization range
         N_points_short=pos_SucClose-Npt_cell;                                % number of discretization points in suction, it is different from N_points       
         distance_out=(theta_vane(end)-theta_SucClose)/(N_points-1);          % angle between each streamtubes in discharge and compression phase
         distance_in=(theta_SucClose)/(N_points-1);                           % angle between each streamtubes in suction phase
         theta_int=0:((theta_SucClose)/(N_points_short-1)):(theta_SucClose);  % discretization of suction phase according to N_points_short
         theta_2=(theta_SucClose):distance_out:theta_vane(end);               % discretization of discharge and compression phase according to number of streamtubes
         theta_1=0:distance_in:(theta_SucClose);                              % discretization of suction phase according to number of streamtubes
         L=2*r*sin((distance_out+distance_in)./4);                            % width of a streamtubes, takes as the mean between inlet and outlet
         
        case 2   %expansion process
            
         N_points=length(theta_vane)-(pos_DisOpen-2*Npt_cell);                          % number of streamtubes, taken as number of points between suction close and end of discretization range
         N_points_short=pos_DisOpen-2*Npt_cell;                                         % number of discretization points in suction, it is different from N_points
         distance_out=(theta_vane(end)-(theta_DisOpen-Gamma))/(N_points-1);             % angular distance between each streamtubes in discharge phase
         distance_in=(theta_DisOpen-Gamma)/(N_points-1);                                % angular distance between each streamtubes in suction and expansion phase
         theta_int=0:((theta_DisOpen-Gamma)/(N_points_short-1)):(theta_DisOpen-Gamma);  % discretization of discharge phase according to number of streamtubes
         theta_2=(theta_DisOpen-Gamma):distance_out:theta_vane(end);                    % discretization of discharge phase according to N_points_short
         theta_1=0:distance_in:(theta_DisOpen-Gamma);                                   % discretization of suction phase according to number of streamtubes
         L=2*r*sin((distance_out+distance_in)./4);                                      % width of a streamtubes, take as the mean between inlet and outlet
    end
   
    delta_theta=theta_2(end:-1:1)-theta_1;
    fact_delta=delta_theta>180;
    delta_theta(fact_delta)=2*pi-delta_theta(fact_delta);
    W_PC=2.*r.*sin(delta_theta./2); % length of each streamtube
    
    %check if distance is lower than vane thickness (s)
    
    flag_W=W_PC<s;
    W_PC(flag_W)=s;
    
    %boundary condition

          f_effective = NaN(size(theta_vane));
    
          fact1 = theta_vane <= theta_SucClose;
          fact2 = theta_vane > theta_SucClose & theta_vane <= (theta_DisOpen-Gamma);
          fact3 = theta_vane > (theta_DisOpen-Gamma) & theta_vane <= (theta_DisClose-Gamma);
          fact4 = theta_vane > (theta_DisClose-Gamma) & theta_vane <= ARC_n_1;
          fact5= theta_vane> ARC_n_1;
          
          switch process
              
              case 1  %compression process
                  
          f_effective(fact1)=0;
          f_effective(fact2)=(((p_comp(length(p_comp)-sum(fact2)+1:end)))-p_in)./(p_out);   
          f_effective(fact3)=1;
          f_effective(fact4) = 1;
          f_effective(fact5) = 1;
          
          P_PC=(f_effective((sum(fact1))+1:end))*(p_out-p_in); % trend of relative pressure in compression and discharge phases
         
          rho_b_pc=mean(rho_b);
          mu_b_pc=mean(mu_b);
          
          Deq= PEclr.*L./(2.*(PEclr + L));    % hydraulic diameter of clearance
          Uo=(P_PC.*0.1./rho_b_pc.*1./(0.35*(4*Deq.*rho_b_pc./mu_b_pc).^(-0.25))./W_PC.*(8.*Deq)).^(1./1.75);  % mean leakage flow velocity
          FlowLeak_PC= Uo.*PEclr.*L.*rho_b_pc;  % leaking mass flow rate in compression and discharge phases

          FlowLeak_PC_resampled=interp1(theta_1,FlowLeak_PC,theta_int,'linear');  % leaking flow rate resampled ofr number of discretization points in suction phases
          FlowLeak_PC_resampled=FlowLeak_PC_resampled.*(sum(FlowLeak_PC)/sum(FlowLeak_PC_resampled));
          
         
          FlowLeak_1=-FlowLeak_PC_resampled(end:-1:1);       % leaking mass flow rate along suction phases (punctual)
          FlowLeak_0=-FlowLeak_PC_resampled(Npt_cell:-1:1);    % leaking mass flow rate along last cell of suction phases  (punctual)
          FlowLeak_2=[FlowLeak_PC, FlowLeak_1(1:Npt_cell)]; % leaking mass flow rate along compression and discharge phases, considering also that part of the last cell will be divided between suction and discharge (punctual)
          
          FlowLeak_suction = -movsum(FlowLeak_1,[0 Npt_cell],'Endpoints','discard');   % leaking mass flow rate along suction phases
          FlowLeak_suction_end=-movsum(FlowLeak_0,[0 Npt_cell],'Endpoints','shrink');      % leaking mass flow rate along last cell of suction phases 
          FlowLeak_comp_dis = -movsum(FlowLeak_2,[0 Npt_cell],'Endpoints','discard');   % leaking mass flow rate along compression and discharge phases, considering also that part of the last cell will be divided between suction and discharge
          FlowLeak=[FlowLeak_suction,FlowLeak_suction_end,FlowLeak_comp_dis];                                % leaking mass flow rate along whole circumference
          
          % Corerection to account for the numerical error linked to the usage of Fourier series instead of the analytical solution
          eps_corr = sum(FlowLeak)/length(FlowLeak);
          FlowLeak = FlowLeak - eps_corr;
          
          %calculation of leakage must be repeated in order to calculate
          %leakage in last cell considering the reduction of volume, so
          %only the contribute of discharge part and not the part of the
          %cell that will be in suction
          FlowLeak_zero = [ FlowLeak_suction,FlowLeak_suction_end, -movsum(FlowLeak_PC,[0 Npt_cell],'Endpoints','shrink')];

              case 2   %expansion process
                  
          f_effective(fact1)=1;
          f_effective(fact2)=(((p_comp(length(p_comp)-sum(fact2)+1:end)))-p_out)./(p_in);   
          f_effective(fact3)=0;
          f_effective(fact4) =0;
          f_effective(fact5) = 0;
          
          p_high=(f_effective(1:(sum(fact2)+sum(fact1))))*((p_in)-p_out);       % trend of relative pressure in suction and expansion phase
          P_PC=interp1(theta_int,p_high,theta_1,'linear');                 % pressure resampled according to number of discretization points in suction
          
          rho_b_pc=mean(rho_b);
          mu_b_pc=mean(mu_b);
          
          Deq= PEclr.*L./(2.*(PEclr + L));    % hydraulic diameter of clearance
          Uo=(P_PC.*0.981./rho_b_pc.*1./(0.35*(4*Deq.*rho_b_pc./mu_b_pc).^(-0.25))./W_PC.*(8.*Deq)).^(1./1.75); % mean leakage flow velocity
          FlowLeak_PC= Uo.*PEclr.*L.*rho_b_pc;
          
          FlowLeak_PC_resampled=interp1(theta_1,FlowLeak_PC,theta_int,'linear');
          FlowLeak_PC_resampled=FlowLeak_PC_resampled.*(sum(FlowLeak_PC)/sum(FlowLeak_PC_resampled));
  
          FlowLeak_1=FlowLeak_PC;       % leaking mass flow rate along suction and expansion phases (punctual)
          FlowLeak_0=FlowLeak_PC(end-Npt_cell:1:end);    % leaking mass flow rate along last cell of expansion phases  (punctual)
          FlowLeak_2=[-FlowLeak_PC_resampled(end:-1:1), FlowLeak_1(1:Npt_cell)]; % leaking mass flow rate along c discharge phases, considering also that part of the last cell will be divided between suction and discharge (punctual)
          
          FlowLeak_suc_exp = -movsum(FlowLeak_1,[0 Npt_cell],'Endpoints','discard');   % leaking mass flow rate along suction and expansion phases
          FlowLeak_exp_end=-movsum(FlowLeak_0,[0 Npt_cell],'Endpoints','shrink');      % leaking mass flow rate along last cell of expansion phases  
          FlowLeak_discharge = -movsum(FlowLeak_2,[0 Npt_cell],'Endpoints','discard');   % leaking mass flow rate along c discharge phases, considering also that part of the last cell will be divided between suction and discharge
          FlowLeak=[FlowLeak_suc_exp,FlowLeak_exp_end(2:end),FlowLeak_discharge];                                % leaking mass flow rate along whole circumference
          
          % Corerection to account for the numerical error linked to the usage of Fourier series instead of the analytical solution
          eps_corr = sum(FlowLeak)/length(FlowLeak);
          FlowLeak = FlowLeak - eps_corr;
          
          %calculation of leakage must be repeated in order to calculate
          %leakage in last cell considering the reduction of volume, so
          %only the contribute of discharge part and not the part of the
          %cell that will be in suction
          FlowLeak_zero = [ FlowLeak_suc_exp,FlowLeak_exp_end(2:end), movsum(FlowLeak_PC_resampled(end:-1:1),[0 Npt_cell],'Endpoints','shrink')];
          
          end

end


end
