function  S7_PrintIO(fileID,IO,FLAG,NUMERIC,PROCESS,GEOMETRY,GAS,LIQUID,LEAK,VANE,SHAFT,BUSHING,NOZZLES,th_flp,th_Vmax,V_inj_nzl,rho_o,nu_o,sigma_o,Dg,fg,epsm,cls_indx,p_out,p_geom,T_g,T_lmean,T_mix,n_pltr,Q_th,L_sp_i,L_sp_m,E_sp_m,eta_g,eta_gT,eta_gl,eta_vol,eta_mecc,m_gas_Id,m_gas,q_gas_Id,q_gas,m_liq_Id,m_liq,Pind_PV,C_shaft_mean,P_mecc,M,Q_g_RS_suction,Q_g_PE_suction,Q_g_VS_suction,Q_g_PE_closed_cell,Q_g_VS_closed_cell,Fleak_g,Fleak_g_PE_suc,Fleak_g_VS_suc,Fleak_g_PE_cell_closed,Fleak_g_VS_cell_closed,Q_g_RS_dis,Q_g_PE_dis,Q_g_VS_dis)
% This function shows on screen the user's inputs and the main results of the simulation 
% INPUT
% fileID [-]      : File identifier of the opened file
% IO [-]          : Structure of basic input output parameters
% FLAG [-]        : structure of flag parameters
% NUMERIC [-]     : Structure of numeric parameters
% PROCESS [-]     : Structure of process parameters
% GEOMETRY [-]    : Structure of geometric parameters
% GAS [-]         : Structure of gas parameters
% LIQUID [-]      : Structure of oil parameters
% LEAK [-]        : Structure of leakage parameters
% VANE [-]        : Structure of vane parameters
% SHAFT [-]       : Structure of shaft parameters
% BUSHING [-]     : Structure of bushing parameters
% NOZZLES [-]     : Structure of nozzles parameters
% th_Vmax [rad]   : port angle that maximize cell volume
% V_inj_nzl [m3]  : volume of oil injected by each nozzle
% rho_o [kg/m3]   : density of injected oil for each nozzles
% nu_o [cSt]      : cinematic viscosity of injected oil for each nozzles
% sigma_o [N/m]   : surface tension of injected oil for each nozzles
% Dg [m]          : droplet diameters for each class
% fg [-]          : droplet frequency for each class
% epsm [-]        : ratio between mass of injected oil and gas mass trapped in a cell for each injector
% cls_indx [-]    : index of disclosed nozzles classes
% p_out [Pa]      : pressure at end of pre-compression 
% p_geom [pa]     : geometrical pressure (pressure in the cell before delivery opening)
% T_g [K]         : gas temperature during closed chamber phase
% T_lmean [K]     : oil adiabatic mixing temperature for each discretization step
% T_mix [K]       : adiabatic mixing temperature during closed chamber phase
% n_pltr [-]      : polytropic index
% Q_th [W]        : Heat power between gas and liquid
% L_sp_i [J/kg]   : specific indicated massic work ( relative to indicated power)
% L_sp_m [J/kg]   : specific mechanical massic work ( relative to mechanical power)
% E_sp_m [J/m3]   : specific mechanical volumetric work 
% eta_g  [%]      : thermodynamic efficiency with respect to isentropic gas compression 
% eta_gT [%]      : thermodynamic efficiency with respect to isothermal compression
% eta_gl [%]      : thermodynamic efficiency with respect to isentropic gas-liquid compression 
% eta_vol [%]     : volumetric efficiency
% eta_mecc [%]    : mechanic efficiency
% m_gas_Id [kg/s] : ideal mass flow rate (no leakages)
% m_gas [kg/s]    : mass flow rate of gas
% q_gas_Id [m3/s] : real gas volumetric flow rate
% q_gas [m3/s]    : ideal gas volumetric flow rate
% m_liq_Id [kg/s] : ideal liquid mass flow rate
% m_liq [kg/s]    : real liquid mass flow rate
% Pind_PV [W] : indicated power
% P_mecc [W]      : mean shaft mechanical power 
% M [-]           : number of injectors

    %% DEFINITIONS %%
    % Titles of report sections
    date   = datestr(IO.Stime,'dd/mm/yyyy HH:MM:SS');
    user   = getenv('username');
    Source = [IO.fullloadname,'.mat'];
    title0 = ['===== SIMULATION REPORT OF ' date ' ======='];
    title1 = 'NUMERIC PARAMETERS';
    title2 = 'PROCESS PARAMETERS';
    title3 = 'GEOMETRY PARAMETERS';
    title4 = 'GAS PARAMETERS'; 
    title5 = 'LIQUID PARAMETERS';
    title6 = 'LEAKAGE PARAMETERS';
    title7 = 'VANE PARAMETERS';
    title8 = 'SHAFT PARAMETERS';
    title9 = 'BUSHING PARAMETERS';
    title10 = 'NOZZLES PARAMETERS & RESULTS';
    title11 = 'RESULTS: THERMODYNAMICS';
    title112 = 'RESULTS: LEAKAGES';
    title12 = 'RESULTS: POWER';
    title13 = 'RESULTS: EFFICIENCY';

    % text dividers
    sep1 = '======================================================';
    sep2 = '------------------------------------------------------';

    % Output formats
    f0  = '%s \n';
    f1  = '%-32s %-10s %-15s \n';
    f2  = '%-32s %-10s %-15.1f \t \n';
    f3  = '%-32s %-10s %-15u \n';
    f4  = '%-32s %-10s %-15.2f \n';
    f5  = '%-32s %-10s %-15.3f \n';
    f6  = '%-32s %-10s ';
    f7  = '%-.0f           \t';
    f8  = '%-.1f       \t';
    f9  = '%-.3f       \t';
    f10 = '%s          \t';
    f11 = '%s       \t';
    f12 = '\n';
    f13 = '%-.1f          \t';
    f14 = '%s            \t';
    f15 = '%-32s %-10s %-15.1e \n';
    f16 = '%-.2f          \t';
    f17 = '%-32s %-10s %-15.5f \n';
    f18 = '%-32s %-9s %-15.5f \n';
    f19 = '%-32s %-5s %-15.5f \n';
    f20  = '%-32s %-9s %-15.3f \n';
    f21  = '%-32s %-5s %-15.3f \n';
    % output preprocessing 
    if PROCESS.process == 1
        process = 'compression';
        used_type = 'compressor'; 
    else
        process = 'expansion';
        used_type = 'expander';
    end
     
    if PROCESS.rot_dir == 1
        rot_dir = 'clockwise'; 
    else
        rot_dir = 'counterclockwise';
    end

    if GEOMETRY.StdProcess == 1
        StdType = 'compressor';
    else
        StdType = 'expander';
    end

    %% PLOT %%
    % Plot report title
    fprintf(fileID,'%s \n',sep1);
    fprintf(fileID,f0,title0);
    fprintf(fileID,'%s \n',sep1);

    % Plot user name
    fprintf(fileID,'USER: %s \n',user);
    fprintf(fileID,'SOURCE: %s \n',Source);

    % Numeric parameters
    fprintf(fileID,'\n%s \n',sep2);
    fprintf(fileID,f0,title1);
    fprintf(fileID,'%s \n',sep2);
    fprintf(fileID,f3, 'Discretization points','-',NUMERIC.Npt_i);
    fprintf(fileID,f2, 'Max iterations','-',NUMERIC.IterMax);
    fprintf(fileID,f2, 'Under-relaxation factor','-',NUMERIC.UndRlx);
    fprintf(fileID,f2, 'Vane discretization step','mm',NUMERIC.stp*1e3);
    fprintf(fileID,f15,'Theorical zero','-',NUMERIC.toll_d);
    fprintf(fileID,f15,'Numerical zero','-',NUMERIC.toll_t);
    % =====================================================================
    % =====================================================================
    fprintf(fileID,f3, 'Shaft discretization points','-',NUMERIC.shm);
    % =====================================================================
    % =====================================================================
    

    % Plot process parameters
    fprintf(fileID,'\n%s \n',sep2);
    fprintf(fileID,f0,title2);
    fprintf(fileID,'%s \n',sep2);
    fprintf(fileID,f1,'Process','-',process);
    fprintf(fileID,f3,'Rotational speed','rpm',PROCESS.rpm);
    fprintf(fileID,f4,'Suction pressure','barA',PROCESS.p_suc*1e-5);
    fprintf(fileID,f4,'Delivery pressure','barA',PROCESS.p_del*1e-5);
    fprintf(fileID,f2,'Inlet gas temperature','°C',PROCESS.T_suc-273.15);
    fprintf(fileID,f2,'Reference temperature','°C',PROCESS.T_0-273.15);
    fprintf(fileID,f4,'Reference pressure','barA',PROCESS.p_0*1e-5);
    fprintf(fileID,f5,'Bushing friction coeff.','-',PROCESS.f_b);
    % =====================================================================
    % =====================================================================
    fprintf(fileID,f1,'Rotation direction','-',rot_dir);
    fprintf(fileID,f2,'angles of x axis w/r/t X','°',rad2deg(PROCESS.zeta));
    % =====================================================================
    % =====================================================================

    % Geometry parameters
    fprintf(fileID,'\n%s \n',sep2);
    fprintf(fileID,f0,title3);
    fprintf(fileID,'%s \n',sep2);
    fprintf(fileID,f1,'Machine','-',GEOMETRY.machine_name);
    fprintf(fileID,f1,'Type','-',StdType);
    fprintf(fileID,f2,'Stator diameter','mm',GEOMETRY.D*1e3); 
    fprintf(fileID,f2,'Rotor diameter','mm',GEOMETRY.d*1e3);
    fprintf(fileID,f2,'Rotor length','mm',GEOMETRY.L*1e3);
    fprintf(fileID,f2,'Vane height','mm',GEOMETRY.l*1e3);
    fprintf(fileID,f4,'Vane thickness','mm',GEOMETRY.s*1e3);
    fprintf(fileID,f3,'Number of vanes','-',GEOMETRY.n_van); 
    fprintf(fileID,f2,'Vane tilt angle','°',rad2deg(GEOMETRY.th_tilt)); 
    fprintf(fileID,f2,'Vane tip radius','mm',GEOMETRY.r_tip*1e3);
    fprintf(fileID,f2,'vane radius offset','mm',GEOMETRY.b*1e3); 
    fprintf(fileID,f2,'Bush diameter','mm',GEOMETRY.d_hub*1e3);
    fprintf(fileID,f2,'Suction open','°',rad2deg(GEOMETRY.th_SucOpen));
    fprintf(fileID,f2,'Suction close','°',rad2deg(GEOMETRY.th_SucClose));
    fprintf(fileID,f2,'Delivery open','°',rad2deg(GEOMETRY.th_DisOpen));
    fprintf(fileID,f2,'Delivery close','°',rad2deg(GEOMETRY.th_DisClose));
    fprintf(fileID,f2,'Rotor-stator clearance','micron',GEOMETRY.RSclr*1e6);
    fprintf(fileID,f2,'Vane-side clearance','micron',GEOMETRY.VSclr*1e6);
    fprintf(fileID,f2,'Rotor-plate clearance','micron',GEOMETRY.PEclr*1e6);
    
    portsCorrection (fileID,PROCESS,GEOMETRY,used_type,StdType,th_flp,th_Vmax,f2)

    % Gas parameters
    fprintf(fileID,'\n%s \n',sep2);
    fprintf(fileID,f0,title4);
    fprintf(fileID,'%s \n',sep2);
    fprintf(fileID,f1,'Gas used','-',GAS.gas_name);
    fprintf(fileID,f2,'specific heat at costant volume','J/kg K',GAS.c_v); 
    fprintf(fileID,f5,'Thermal conductivity','W/m K',GAS.k_g);
    fprintf(fileID,f5,'Molecolar mass','kg/kmol',GAS.MM_g*1e3);
    fprintf(fileID,f6,'Gas composition','-');
    for i = 1:size(GAS.molcomp_g,2)
    fprintf(fileID,f11,GAS.name4prop{i});       
    end
    fprintf(fileID,f12);
    fprintf(fileID,f6,'Molar fraction','-');
    for i = 1:size(GAS.molcomp_g,2)
    fprintf(fileID,f9,GAS.molcomp_g(i));       
    end
    fprintf(fileID,f12);        
    fprintf(fileID,f1,'Thermodynamic model','-',GAS.model_g);        
    fprintf(fileID,f1,'Gas property source','-',GAS.propSource);  


    % Oil parameters
    fprintf(fileID,'\n%s \n',sep2);
    fprintf(fileID,f0,title5);
    fprintf(fileID,'%s \n',sep2);
    fprintf(fileID,f1,'Oil used','-',LIQUID.oil_name);
    fprintf(fileID,f2,'Reference temperature','°C',LIQUID.T_ref-273.15);
    fprintf(fileID,f2,'Density at reference temperature','kg/m3',LIQUID.rho_ref);
    fprintf(fileID,f2,'Specific heat','J/kg K',LIQUID.c_l);        
    fprintf(fileID,f5,'Thermal conductivity','W/m K',LIQUID.k_l); 
    fprintf(fileID,f2,'Kinemaatic viscosity at 40°C','cSt',LIQUID.nu_40); 
    fprintf(fileID,f2,'Kinematic viscosity at 100°C','cSt',LIQUID.nu_100); 
    
    % Plot leakages parameters
    fprintf(fileID,'\n%s \n',sep2);
    fprintf(fileID,f0,title6);
    if FLAG.fLKG == 1
    fprintf(fileID,f1,'Leakage analysys','-','Enabled');    
    fprintf(fileID,f1,'Leakage model rotor-stator (RS)','-',LEAK.RSmodel);
    fprintf(fileID,f1,'Two-phase viscosity model RS','-',LEAK.MUmodel_RS);
    fprintf(fileID,f1,'Leakage model vane-side (VS)','-',LEAK.VSmodel);
    fprintf(fileID,f1,'Two-phase viscosity model VS','-',LEAK.MUmodel_VS);
    fprintf(fileID,f1,'Leakage model rotor-plate (PE)','-',LEAK.PEmodel);
    fprintf(fileID,f1,'Two-phase viscosity model PE','-',LEAK.MUmodel_PE);
    else
    fprintf(fileID,f1,'Leakage analysys','-','Disabled');    
    end

    % Plot vane parameters
    fprintf(fileID,'\n%s \n',sep2);
    fprintf(fileID,f0,title7);
    fprintf(fileID,'%s \n',sep2);
    fprintf(fileID,f1,'Vane type','-',VANE.vane_type);
    fprintf(fileID,f4,'Vane density','kg/m3',VANE.rho_vane);
    fprintf(fileID,f2,'Vane G offset','%',VANE.offsetG); 
    fprintf(fileID,f5,'Slot friction coeff.','-',VANE.f_c);
    fprintf(fileID,f5,'Tip friction coeff.','-',VANE.f_t);
    
    % =====================================================================
    % =====================================================================
    % Plot shaft parameters
    fprintf(fileID,'\n%s \n',sep2);
    fprintf(fileID,f0,title8);
    fprintf(fileID,'%s \n',sep2);
    fprintf(fileID,f4,'Rotor+shaft density','kg/m3',SHAFT.rho_s);
    
    % Plot bushing parameters
    fprintf(fileID,'\n%s \n',sep2);
    fprintf(fileID,f0,title9);
    fprintf(fileID,'%s \n',sep2);
    fprintf(fileID,f2,'Rotor end-bushing distance','mm',BUSHING.d1*1e3);
    fprintf(fileID,f2,'Bushing length','mm',BUSHING.bB*1e3);
    % =====================================================================
    % =====================================================================

    % Plot Nozzles parameters
    fprintf(fileID,'\n%s \n',sep2);
    fprintf(fileID,f0,title10);
    fprintf(fileID,'%s \n',sep2);
    fprintf(fileID,f3,'Nusselt Number','-',NOZZLES.Nu);
    fprintf(fileID,f6,'Nozzle type','-');
    for i=1:M
        fprintf(fileID,f10,NOZZLES.type_nz{i});
    end
    fprintf(fileID,f12);

    fprintf(fileID,f6,'Nozzles name','-');
    for i=1:M
        fprintf(fileID,f10,NOZZLES.name_nz{i});
    end
    fprintf(fileID,f12);
    fprintf(fileID,f6,'Nozzle position','°');
    fprintf(fileID,f8,rad2deg(NOZZLES.theta_nz)');
    fprintf(fileID,f12);
    fprintf(fileID,f6,'Nozzle number','-');
    fprintf(fileID,f7,NOZZLES.num_nz');
    fprintf(fileID,f12);
    fprintf(fileID,f6,'Oil injection pressure','barA');
    fprintf(fileID,f9,NOZZLES.p_inj*10^-5);
    fprintf(fileID,f12);
    fprintf(fileID,f2,'Inlet oil temperature','°C');
    fprintf(fileID,f13,NOZZLES.Tl_in-273.15);
    fprintf(fileID,f12);
    fprintf(fileID,f6,'Nozzle number of classes','-');
    fprintf(fileID,f7,NOZZLES.num_cls');
    fprintf(fileID,f12);
    fprintf(fileID,f6,'Excluded droplet distribution','-');
    fprintf(fileID,f16,NOZZLES.alpha_nz');
    fprintf(fileID,f12);

    fprintf(fileID,f6,'Orefice diameter','mm');
    for i=1:M
        if isfield(NOZZLES.specs{i},'D_or')
            fprintf(fileID,f13, NOZZLES.specs{i}.D_or);
        else
            fprintf(fileID,f14,'-');
        end
    end
    fprintf(fileID,f12);

    fprintf(fileID,f6,'Orefice length','mm');
    for i=1:M
        if isfield(NOZZLES.specs{i},'L_or')
            fprintf(fileID,f13, NOZZLES.specs{i}.L_or);
        else
            fprintf(fileID,f14,'-');
        end
    end
    fprintf(fileID,f12);

    fprintf(fileID,f6,'Reference flow rate','l/min');
    for i=1:M
        if isfield(NOZZLES.specs{i},'V_ref')
            fprintf(fileID,f13, NOZZLES.specs{i}.V_ref);
        else
            fprintf(fileID,f14,'-');
        end
    end
    fprintf(fileID,f12);

    fprintf(fileID,f6,'Reference pressure drop','bar');
    for i=1:M
        if isfield(NOZZLES.specs{i},'Delta_p_ref')
            fprintf(fileID,f13, NOZZLES.specs{i}.Delta_p_ref/10^5);
        else
            fprintf(fileID,f14,'-');
        end
    end
    fprintf(fileID,f12);

    fprintf(fileID,f6,'Spray angular aperture','°');
    for i=1:M
        if isfield(NOZZLES.specs{i},'gamma_nz')
            fprintf(fileID,f13, NOZZLES.specs{i}.gamma_nz);
        else
            fprintf(fileID,f14,'-');
        end
    end
    fprintf(fileID,f12);
    fprintf(fileID,f6,'Oil injected Volume','mm3');
    fprintf(fileID,f9,V_inj_nzl(1:end)*10^9);
    fprintf(fileID,f12);
    fprintf(fileID,f6,'Oil/gas mass ratio','-');
    fprintf(fileID,f9,epsm);
    fprintf(fileID,f12);
    fprintf(fileID,f6,'Oil density','kg/m3');
    fprintf(fileID,f9,rho_o(1:end));
    fprintf(fileID,f12);
    fprintf(fileID,f6,'Oil viscosity','cSt');
    fprintf(fileID,f9,nu_o(1:end));
    fprintf(fileID,f12);
    fprintf(fileID,f6,'Surface tension','N/m');
    fprintf(fileID,f9,sigma_o(1:end));
    fprintf(fileID,f12);

    fprintf(fileID,f6,'Droplets diameter per class','µm');
    fprintf(fileID,f7,(Dg(cls_indx==1)*10^6)');
    fprintf(fileID,f12);
    for i=2:M
        fprintf(fileID,f6,'','');
        fprintf(fileID,f7,(Dg(cls_indx==i)*10^6)');
        fprintf(fileID,f12);
    end

    fprintf(fileID,f6,'Class frequency','-');
    fprintf(fileID,f9,(fg(cls_indx==1))');
    fprintf(fileID,f12);
    for i=2:M
        fprintf(fileID,f6,'','');
        fprintf(fileID,f9,(fg(cls_indx==i))');
        fprintf(fileID,f12);
    end
    
    % Plot thermoynamics results
    fprintf(fileID,'\n%s \n',sep2);
    fprintf(fileID,f0,title11);
    fprintf(fileID,'%s \n',sep2);
    fprintf(fileID,f4,'Delivery pressure','bar',p_out*1e-5);
    fprintf(fileID,f4,'Delivery pressure (geometrical)','bar',p_geom*1e-5);
    fprintf(fileID,f5,'Gas mass flow rate (ideal)','kg/s',m_gas_Id);
    fprintf(fileID,f5,'Gas mass flow rate (actual)','kg/s',m_gas);
    fprintf(fileID,f18,'Gas volum. flow rate suction (RS)','kg/s',Fleak_g);
    fprintf(fileID,f18,'Gas volum. flow rate suction (PE)','kg/s',Fleak_g_PE_suc);
    fprintf(fileID,f18,'Gas volum. flow rate suction (VS)','kg/s',Fleak_g_VS_suc);
    fprintf(fileID,f19,'Gas volum. flow rate closed cell (PE)','kg/s',Fleak_g_VS_cell_closed);
    fprintf(fileID,f19,'Gas volum. flow rate closed cell (VS)','kg/s',Fleak_g_VS_cell_closed);
    fprintf(fileID,f5,'Gas volum. flow rate (ideal)','l/min',q_gas_Id*1e3*60);
    fprintf(fileID,f5,'Gas volum. flow rate (actual)','l/min',q_gas*1e3*60);
    fprintf(fileID,f5,'Gas volum. flow rate (RS)','l/min',Q_g_RS_dis*1e3*60);
    fprintf(fileID,f5,'Gas volum. flow rate (PE)','l/min',Q_g_PE_dis*1e3*60);
    fprintf(fileID,f5,'Gas volum. flow rate (VS)','l/min',Q_g_VS_dis*1e3*60);
    fprintf(fileID,f5,'Oil mass flow rate (ideal)','kg/s',m_liq_Id);
    fprintf(fileID,f5,'Oil mass flow rate (actual)','kg/s',m_liq);
    fprintf(fileID,f2,'Outlet gas temperature','°C',T_g(end)-273.15);
    fprintf(fileID,f2,'Mean outlet oil temperature','°C',T_lmean(end)-273.15);
    fprintf(fileID,f2,'Outlet adiabatic mix temperature','°C',T_mix(end)-273.15);
    fprintf(fileID,f5,'Heat exchanged gas/oil','kW',Q_th*1e-3);
    
    if FLAG.fLKG_plot
    % Plot leakages results
    fprintf(fileID,'\n%s \n',sep2);
    fprintf(fileID,f0,title112);
    fprintf(fileID,'%s \n',sep2);
    fprintf(fileID,f20,'Gas volum. flow rate suction (RS)','l/min',Q_g_RS_suction*1e3*60);
    fprintf(fileID,f20,'Gas volum. flow rate suction (PE)','l/min',Q_g_PE_suction*1e3*60);
    fprintf(fileID,f20,'Gas volum. flow rate suction (VS)','l/min',Q_g_VS_suction*1e3*60);
    fprintf(fileID,f21,'Gas volum. flow rate closed cell (VS)','l/min',Q_g_VS_closed_cell*1e3*60);
    fprintf(fileID,f21,'Gas volum. flow rate closed cell (PE)','l/min',Q_g_PE_closed_cell*1e3*60);
    end
    
    % Plot power results
    fprintf(fileID,'\n%s \n',sep2);
    fprintf(fileID,f0,title12);
    fprintf(fileID,'%s \n',sep2);
    fprintf(fileID,f5,'Indicated power','kW',Pind_PV*1e-3);
    fprintf(fileID,f5,'Mechanical torque','Nm',C_shaft_mean);
    fprintf(fileID,f5,'Mechanical power','kW',P_mecc*1e-3);
    fprintf(fileID,f4,'Spefific indicted work','kJ/kg',L_sp_i*1e-3);
    fprintf(fileID,f4,'Specific mechanical work','kJ/kg',L_sp_m*1e-3);
    fprintf(fileID,f5,'Specific energy (by mech power)','kW*min/m3',E_sp_m*1e-3/60);

    % Plot efficiency results
    fprintf(fileID,'\n%s \n',sep2);
    fprintf(fileID,f0,title13);
    fprintf(fileID,'%s \n',sep2);
    fprintf(fileID,f4,'Isoentropic efficiency','%',eta_gl); 
    fprintf(fileID,f4,'Mechanic efficiency','%',eta_mecc);
    fprintf(fileID,f4,'Gas-only process efficiency','%',eta_g); 
    fprintf(fileID,f4,'Efficiency respect to isotherm','%',eta_gT);
    fprintf(fileID,f4,'Volumetric efficiency','%',eta_vol);
    %exclude NaN values from politropic index array
    goodn_pltr=~isnan(n_pltr);
    fprintf(fileID,f4,'Mean politropic index','-',mean(n_pltr(goodn_pltr))); 
    fprintf(fileID,'\n%s \n',sep2);
    
    % if stress analysis is not performed, report ends here; computational
    % time is plotted
    if FLAG.fSTR == 0
    fprintf(fileID, 'Elapsed time: %.2f s \n',toc(IO.CALLtime));
    end
end

function portsCorrection (fileID,PROCESS,GEOMETRY,used_type,StdProcess,th_flp,th_Vmax,f2)
% This function display the change in ports angle due to flipping or volume maximization
%  INPUT
% fileID [-]         : File identifier of the opened file
% PROCESS [-]        : Structure of process parameters
% GEOMETRY           : Structure of geometric parameters
% used_type [string] : Process performed    
% StdProcess         : Design process
% th_flp             : Array of flipped ports angle
% th_Vmax            : Angle changed to maximize volume
% f2                 : fprintf options

    fFLIP = 0;
    fMAX  = 0;
    sep = '******************';
    
    if PROCESS.process ~= GEOMETRY.StdProcess
        fFLIP = 1;
    end
    if isnan(GEOMETRY.th_SucClose) || isnan(GEOMETRY.th_DisOpen)
        fMAX = 1;
    end
    
    if fFLIP || fMAX
        fprintf(fileID,'\n%s ',sep);
        fprintf(fileID,'PORTS CORRECTION');
        fprintf(fileID,' %s \n',sep);
    end

    if fFLIP
        fprintf(fileID,'The machine is: %s, but is used as: %s.\n', upper(StdProcess),upper(used_type));
        fprintf(fileID,'Ports angles must be flipped.\n');
        fprintf(fileID,'New ports angles:\n\n');
        fprintf(fileID,f2,'Suction open','°',rad2deg(th_flp(1)));
        fprintf(fileID,f2,'Suction close','°',rad2deg(th_flp(2)));
        fprintf(fileID,f2,'Delivery open','°',rad2deg(th_flp(3)));
        fprintf(fileID,f2,'Delivery close','°',rad2deg(th_flp(4)));
    end

    if fMAX
        fprintf(fileID,'Volume maximization has been performed:\n');
        switch PROCESS.process
            case 1
                fprintf(fileID,'Suction close angle maximizes suction volume in %s.\n', used_type);
                fprintf(fileID,'New ports angles:\n\n');
                fprintf(fileID,f2,'Suction close','°',rad2deg(th_Vmax));
            case 2
                fprintf(fileID,'Delivery open angle maximizes delivery volume in %s.\n', used_type);
                fprintf(fileID,'New ports angles:\n\n');
                fprintf(fileID,f2,'Delivery open','°',rad2deg(th_Vmax));
        end
    end
    fprintf(fileID,'******************************************************\n');
end