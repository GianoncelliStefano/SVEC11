function CALL_SVEC_M135A
    % Works with Matlab version 2018a or more recent version.
    IO.CALLtime = tic;  % start simulation timer     
    IO.Stime    = now;  % current date
    
    %% INPUT - OUTPUT %%
    % If you want to execute only postprocessing from saved data:
    %  * Set IO.fMODE = 0
    %  * Specify name and path of simulation file to load
    %  * Edit this section with your preferences
    %  * Run this code

    IO.save_dir = 'SVEC Results';                             % Directory    - default : '..\SVEC Results';
    IO.filename = datestr(now,'yyyy_mm_dd  HH_MM_SS');           % Filename     - default : datestr(now,'yyyy_mm_dd  HH_MM_SS');
    IO.fMODE    = 1;      % Save\Load mode     0: Load Simul     1: New Simul

    IO.fDRPT    = 1;      % Display Report     0: Disable        1: Enable
    IO.fSRPT    = 0;      % Save Report        0: Disable        1: Enable
    IO.fDPLT    = 0;      % Display Plots      0: Disable        1: Enable
    IO.fSPLT    = 0;      % Save Plots         0: Disable        1: Enable
    IO.ext      = 'fig';  % Extension to save plot (see 'saveas' function for more)
                          % | 'fig' | pdf | 'eps' (b&w) | 'epsc' (color) | 'png' | 'jpeg' | ...

    [IO,fOK] = SX_Startup (IO);

    if IO.fMODE         % New simulation
        %% FLAG %%
        FLAG.fDBG = 1;      % Controls                  0: User only      1: User + Developer
        FLAG.fSTR = 0;      % Stress Analysis;          0: Disable        1: Enable
        FLAG.fLKG = 1;      % Leakages                  0: Disable        1: Enable
        FLAG.fLKG_in = 1;   % Leakages path 2 and 3     0: Disable        1: Enable
        FLAG.fLKG_plot =0;  % plot Leakages             0: Disable        1: Enable

        %% NUMERIC %%
        NUMERIC.Npt_i   = 10000;   % [-]  Number of grid discetization points     (Reccomended MIN: 10000 MAX: 100000)
                                   % NOTE: in case of activated Stress Analysis
                                   % it is required to set Npt_i = 99000
                                   % in order to achieve a satifying
                                   % resolution together with a proper computational
                                   % time. A resampling method is implemented in such case.                                  
        NUMERIC.IterMax = 30;      % [-]  Maximum number of iteration
        NUMERIC.UndRlx  = 0.9;     % [-]  Under-relaxation factor for iteration   (Reccomended MIN: 0.5 MAX: 1,5) (if low: slow and stable convergnece   if high: fast and instable convergence)
        NUMERIC.stp     = 0.1;     % [mm] Discretization step for stress analysis (Reccomended MIN: 0.1   MAX: 1)
        NUMERIC.toll_d  = 1e-3;    % [-]  Numeric tolerance (default: 1e-3)
        NUMERIC.toll_t  = 1e-9;    % [-]  Theoric tolerance (default: 1e-9)
        
        % ==================================================================
        % ==================================================================
        NUMERIC.shm     = 500;     % [-] Rotor+shaft axial discretization points
        % ==================================================================
        % ==================================================================

        %% PROCESS %%
        % Creation of PROCESS structure with operating conditions
        PROCESS.process = 1;        % [bool]  Process performed    1: compression        2: Expansion
        PROCESS.rpm     = 1500;     % [rpm]   Rotational speed
        PROCESS.p_suc   = 1;     % [barA]  Suction pressure
        PROCESS.p_del   = 8.5;        % [barA]  Delivery pressure - (use NaN for geometric pressure)
        PROCESS.T_suc   = 28;     % [°C]    Gas inlet temperature
        PROCESS.T_0     = 20.0;     % [°C]    Reference temperature for energy balance
        PROCESS.p_0     = 1;       % [barA]  Reference pressure
        PROCESS.f_b     = 0.008;    % [adim]  Bushing friction coefficient
        % ==================================================================
        % ==================================================================
        PROCESS.rot_dir = 'clockwise';       % direction of rotation: 'clockwise' or 'counterclockwise'
        PROCESS.zeta = 0;                    % [°] angles of x axis with respect to X (absolute system of reference). 
                                             % x is always oriented in the opposite
                                             % direction of tangency point.
        % ==================================================================
        % ==================================================================


        %% GEOMETRY %%
        % Creation of GEOMETRY structure: manual input (setting Machine_selected = 'UserDefined' and GEOMETRY.c = 1 or 2)
        % or from database (e.g. setting Machine_selected = 'M111H')
        % NOTE: * The ports angles must always be inserted so that the suction open is the lowest angle and the delivery close is the highest angle,
        %           indipendently from the process chosen (compression or expansion)
        %       * Select the correct design process of the machine
        %       * If GEOMETRY.StdProcess and PROCESS.process are different, SVEC will flip the angles of the machine
        %
        %               AVAILABLE COMPRESSORS                  AVAILABLE EXPANDERS
        % ||  'XT65'  |  'M80F'  | 'M86G'   |  'M195B'   || 'M100Xorc' (A) |           ||
        % ||  'M65C'  |  'M86B'  | 'M86J'   |  'M195C'   || 'M100Xtpg' (B) |           ||
        % ||  'M75'   |  'M86C'  | 'M111H'  |  'M195BX'  || 'M200X'        |           ||
        % ||  'M80B'  |  'M86D'  | 'M135H'  |  'M111C'   ||                |           ||
        % ||  'M80D'  |  'M86E'  | 'M170K'  |  'M135L'   ||                |           ||
        % ||  'M215'  |  'M135B' | 
        Machine_selected = 'M135A';

        if strcmp(Machine_selected,'UserDefined')
            GEOMETRY.c = 1;     % stator geometry     1: circular   2: elliptical stator;

            switch GEOMETRY.c
                case 1   % cylindrical geometry
                    GEOMETRY.StdProcess  = 1;      % [-]  design process for this machine  1: Compression   2: Expansion
                    GEOMETRY.D           = 270;   % [mm] stator diameter
                    GEOMETRY.d           = 226.12; % [mm] rotor diameter
                    GEOMETRY.L           = 452.24; % [mm] rotor length
                    GEOMETRY.l           = 66.5;   % [mm] vane height
                    GEOMETRY.s           = 8;      % [mm] vane thickness (for s<toll_t, 1D vane geometry model is used)
                    GEOMETRY.s_1D        = 3.96;   % [mm] vane thickness for 1D geometry
                    GEOMETRY.n_van       = 7;      % [-]  number of vanes
                    GEOMETRY.d_hub       = 81.17;  % [mm] bushing diameter
                    GEOMETRY.th_SucOpen  = 48.5;   % [°]  suction opening angle
                    GEOMETRY.th_SucClose = 161.6;  % [°]  suction close angle     (COMPRESSION: use NaN to use  maxize suction volume)
                    GEOMETRY.th_DisOpen  = 325.3;  % [°]  delivery opening angle  (EXPANSION:   use NaN to maximize delivery volume)
                    GEOMETRY.th_DisClose = 356.3;  % [°]  delivery closing angle
                    GEOMETRY.th_tilt     = 0;      % [°]  vane tilt angle positive/null/negative for forward/radial/backward vanes
                    GEOMETRY.b           = 0;      % [mm] offset between center of tip circonference and vane axis
                    GEOMETRY.r_tip       = 9;      % [mm] vane tip radius
                    GEOMETRY.TgAngle     = 0;      % [°]  sealing arc between stator and rotor
                    GEOMETRY.RSclr       = 50;     % [micron] clearance between rotor and stator
                case 2  % elliptical geometry
                    GEOMETRY.StdProcess  = 1;      % [-]  design process for this machine  1: Compression   2: Expansion
                    GEOMETRY.d           = 111;    % [mm] rotor diameter
                    GEOMETRY.e           = 0.6;    % [-]  ellipse eccentricity
                    GEOMETRY.L           = 275;    % [mm] rotor length
                    GEOMETRY.l           = 38;     % [mm] vane height
                    GEOMETRY.s           = 4.72;   % [mm] vane thickness
                    GEOMETRY.n_van       = 8;      % [-]  numebr of vanes
                    GEOMETRY.d_hub       = 30;     % [mm] bushing diameter
                    GEOMETRY.th_SucOpen  = 0;      % [°] suction opening angle
                    GEOMETRY.th_SucClose = 65;     % [°] suction close angle     (COMPRESSION: use NaN to maxize suction volume)
                    GEOMETRY.th_DisOpen  = 150;    % [°] delivery opening angle  (EXPANSION:   use NaN to maximize delivery volume)
                    GEOMETRY.th_DisClose = 170;    % [°] delivery closing angle
                    GEOMETRY.TgAngle     = 0;      % [°] sealing arc between stator and rotor (TgAngle = 0 if no socket is present)
                    GEOMETRY.RSclr       = 50;     % [micron] clearance between rotor and stator
                otherwise
                    warning('CALL_SVECmodelMain:Logic','Invalid value of GEOMETRY.c');
                    SX_Logfile ('e',{lastwarn});
                    fOK = 0;
            end
        else
            [GEOMETRY,fOK] = SX_DatabaseLoad('Mach', Machine_selected,fOK);
        end
        GEOMETRY.VSclr = 225;                            % [micron] vane-end-plate clearance size
        GEOMETRY.PEclr = 225;                            % [micron] rotor-end-plate clearance size
        GEOMETRY.machine_name = Machine_selected;
        GEOMETRY.RSclr = 72.5;
        clear Machine_selected
        
        %% GAS %%
        % Creation of GAS structure: manual input (setting gas_selected = 'UserDefined')
        % or from database (e.g. setting gas_selected = 'CO2')
        %
        %                        AVAILABLE GAS
        % ||   'Air'     | 'Hydrogen' |  'Methane'  |  'R236fa'   ||
        % ||   'Argon'   | 'Nitrogen' |  'Steam'    |  'R1233zd'  ||
        % ||   'Helium'  | 'Xenon'    |  'CO2'      |             ||
        Gas_selected   = 'Air';

        if strcmp(Gas_selected,'UserDefined')
            GAS.c_v        = 720;                 % [J/kg K]  specific heat at costant volume
            GAS.k_g        = 0.028;               % [W/ m K]  thermal conductivity
            GAS.MM_g       = 28.959;              % [kg/kmol] molecolar mass
            GAS.name4prop  = {'Methane','Steam'}; % [-]       name of species contained in gas (See PureCompData.m)
            GAS.molcomp_g  = [0.72 , 0,28];       % [-]       molar fraction of contained gas
        else
            [GAS,fOK]      = SX_DatabaseLoad('Gas',Gas_selected,fOK);
        end

        GAS.model_g    = 'simple';    % 'simple' : simplified thermodynamics | 'full': full thermodynamics
        GAS.propSource = 'gas';       % 'gas' : ThermoPhysProps, but with imposed gas phase(NO resolution of VLE) | 'direct' : ThermoPhysProps | 'REFPROP': link to REFPROP
        GAS.model_full = 'fullTP';
        if strcmp(GAS.model_g,'full')
            current_directory = cd;
            tpp_directory     = strcat(current_directory,"\ThermoPhysProps");
            addpath(tpp_directory);
        end
        
        GAS.gas_name   = Gas_selected;
        clear Gas_selected

        %% LIQUID %%
        % Creation of LIQUID structure: manual input (setting Liquid_selected = 'UserDefined')
        % or from database (e.g. setting Liquid_selected = 'F2')
        %
        %             AVAILABLE LIQUID
        % ||   'F2'   |  'AND5010'  |  'CPI4200'  ||
        % ||   'F4'   |  'CPI4201'  |             ||
        Liq_selected = 'F2';

        if strcmp(Liq_selected,'UserDefined')
            LIQUID.T_ref   = 15;    % [°C]     reference temperature
            LIQUID.rho_ref = 1500;  % [kg/m3]  density at reference temperature
            LIQUID.c_l     = 1940;  % [J/kg K] specific heat
            LIQUID.k_l     = 0.128; % [W/ m K] thermal conductivity
            LIQUID.nu_40   = 50;    % [cSt]    kinemaatic viscosity at 40°C
            LIQUID.nu_100  = 11;    % [cSt]    kinematic viscosity at 100°C
        else
            [LIQUID,fOK] = SX_DatabaseLoad('Oil',Liq_selected,fOK);
        end
        LIQUID.oil_name = Liq_selected;
        clear Liq_selected

       %% LEAKAGES %%
        LEAK.MUmodel_RS ='Awad';             % two-phase viscosity model
                                                %  | 'Dukler' (low viscosity)| | 'Maxwell' (High viscosity)| | 'Awad' |
        LEAK.RSmodel    ='Poiselle-Couette';    % Leakage model for rotor stator clearance
                                                %  | 'Poiselle-Couette' | 'Ferreira' | 'Gasche' |
        LEAK.MUmodel_VS ='Awad';             % two-phase viscosity model
                                                %  | 'Dukler' (low viscosity)| | 'Maxwell' (High viscosity)| | 'Awad' |   
        LEAK.VSmodel    ='Poiselle-Couette';    % Leakage model for vane side end-plates clearance
                                                %  | 'Poiselle-Couette' | 'Yuan' |    
        LEAK.MUmodel_PE ='Awad';              % two-phase viscosity model
                                                %  | 'Dukler' (low viscosity)| | 'Maxwell' (High viscosity)| | 'Awad' |
        LEAK.PEmodel    ='Poiselle-Couette';                % Leakage model for plate-ends clearance
                                                %  | 'Badr' |  |' Poiselle-Couette' |
   
        %% VANE %%
        % Creation of VANE structure: manual input (setting vane_selected = 'UserDefined')
        % or from database (e.g. setting Liquid selected = 'standard')
        %
        %                  AVAILABLE VANE
        % ||  'Standard'   |  'Tenmat'  |  'Sandwich'   ||
        %
        % NOTE: offset percentage is referred to half of vane length (l/2).
        %       e.g.: if offsetG = 10% means that center of mass is closer to vane
        %       tip by a value of (l/2)*10/100 respect to geometric center of mass.
        %       As a consequence:
        %       - offsetG =    0%   center of mass correspond to geometric one
        %       - offsetG = + 100%  center of mass is located on vane tip
        %       - offsetG = - 100%  center of mass is located on vane bottom
        %
        % M195B vanes are drilled, so their mass is reduced. Multiply rho_van*0.5416

        Vane_selected = 'Standard';

        if strcmp(Vane_selected,'UserDefined')
            VANE.rho_vane  = 7000;   % [kg/m3] vane density
            VANE.offsetG   = 10;     % [%]     center of mass offset (see NOTE)
            VANE.f_c       = 0.12;   % [-]     slot friction coefficient
            VANE.f_t       = 0.03;   % [-]     tip friction coefficient
        else
            [VANE,fOK] = SX_DatabaseLoad('Vane',Vane_selected,fOK);
        end
        VANE.vane_type = Vane_selected;
        clear Vane_selected

        %% SHAFT %%
        % ================================================================
        % ================================================================
        % Creation of SHAFT structure
        SHAFT.rho_s = 7100;     % [kg/m3] Rotor+shaft material density 
        % ================================================================
        % ================================================================
        %% BUSHING %%
        % ================================================================
        % ================================================================
        % Creation of BUSHING structure
        BUSHING.d1 = 5;                          % [mm] Rotor end-bushing distance 
        BUSHING.bB = 60;                         % [mm] Bushing length 
        % ================================================================
        % ================================================================
        %% NOZZLES %%
        % Nozzles with the lowest angle will be the closest to the suction of the machine,
        % indipendently from the process performed (compression or expansion)
        %
        %                    AVAILABLE NOZZLES
        % || 'RD734' (fc)    | 'RD816' (hc) |  'RD817' (hc) ||
        % || 'ud' (po/hc/fc) | 'un' (un)    |               ||
        
        NOZZLES.f_nz     = [     0;        0;        1;        0;        0];  % [bool]   nozzles status 1: enabled    0: disabled
        NOZZLES.theta_nz = [ 213.2;    233.9;      247;    255.3;    277.4];  % [°]      nozzles angular position
        NOZZLES.num_nz   = [     2;        2;        6;        2;        2];  % [-]      numer of nozzles for each angular position
        NOZZLES.p_inj    = [  6.80;     6.80;       7.5;     6.80;     6.80];  % [barA]   injection pressure
        NOZZLES.Tl_in    = [  65.0;     65.0;     80.0;     65.0;     65.0];  % [°C]    injection temperature
        NOZZLES.num_cls  = [     7;        7;        1;        7;        7];  % [-]      number of classes for each nozzles
        NOZZLES.alpha_nz = [  0.01;     0.01;     0.01;     0.01;     0.01];  % [-]      drop distribution excluded during discretization (must be between 0 and 1 - default value: 0.01)
        NOZZLES.type_nz  = {  'fc';     'fc';     'po';     'fc';     'fc'};  % [string] type of nozzles
                                                                                % 'po' plain orifice nozzle
                                                                                % 'fc' full-cone pressure-swirl nozzle
                                                                                % 'hc' hollow-cone pressure swirl nozzle
                                                                                % 'kf' known flow nozzle
                                                                                % 'un' undefined nozzle
        NOZZLES.name_nz  = {'RD734'; 'RD734';    'ud';    'RD734'; 'RD734'};  % [string] nozzles name
                                                                                % chose a name from MATTEI database or:
                                                                                % 'un' for undefined nozzle
                                                                                % 'ud' user-defined nozzle
        NOZZLES.Nu       = 2;                                    % [-] Nusselt number (default value: 2)

        % For UNDEFINED 'un' nozzles create a structure NOZZLES.specs with the fields:
        %  D_g   [m]  mean drop diameter for each class
        %  f_g   [-]  frequency for each class
        %  eps_m [-]  oli/gas mass ratio for each nozzle
        %  NOTE: For undefined nozzle in suction, eps_m is defined respect to ideal
        %        gas mass in first cloced cell
        % Insert between curly brackets {} the index of the corresponding nozzle
        % e.g.:
        % NOZZLES.specs{1}.D_g    = [300, 400].*1e-6;
        % NOZZLES.specs{1}.f_g    = [0.5, 0.5];
        % NOZZLES.specs{1}.eps_m  = 0.05;

        % For USER DEFINED 'ud' KNOWN FLOW nozzles create a structure NOZZLES.specs with the fields:
        %  D_g   [m]     droplet diameters for each class
        %  f_g   [-]     frequency for each class
        %  V_inj [m3/s] total injected volume
        % Insert between curly brackets {} the index of the corresponding nozzle
        % e.g.:
        % NOZZLES.specs{2}.D_g   = [100, 200, 300, 400].*1e-6;
        % NOZZLES.specs{2}.f_g   = [0.15, 0.25, 0.35, 0.0.25];
        % NOZZLES.specs{2}.V_inj = [90]*60*1000;
        
        % For USER DEFINED 'ud' PLAIN ORIFICE nozzles create a structure NOZZLES.specs with the fields:
        %  D_or [mm]  orifice diameter
        %  L_or [mm]  orifice length
        % Insert between curly brackets {} the index of the corresponding nozzle
        % e.g.:
         NOZZLES.specs{3}.D_or = 2.9;
         NOZZLES.specs{3}.L_or = 14.49;

        % For USER DEFINED 'ud' FULL CONE PRESSURE-SWIRL nozzles create a structure NOZZLES.specs with the fields:
        %  D_or [mm]        final nozzle diameter
        %  V_ref [l/min]    water flow rate supplied with the reference pressure drop
        %  Delta_p_ref [Pa] reference pressure drop
        % Insert between curly brackets {} the index of the corresponding nozzle
        % e.g.:
        % NOZZLES.specs{1}.D_or        = 2.8;
        % NOZZLES.specs{1}.V_ref       = 3.5;
        % NOZZLES.specs{1}.Delta_p_ref = 4*10^5;

        % For USER DEFINED 'ud' HOLLOW CONE PRESSURE-SWIRL nozzles create a structure NOZZLES.specs with the fields:
        %  D_or [mm]        final nozzle diameter
        %  V_ref [l/min]    water flow rate supplied with the reference pressure drop
        %  Delta_p_ref [Pa] reference pressure drop
        %  gamma_nz [°]     angular aperture of cone
        % Insert between curly brackets {} the index of the corresponding nozzle
        % e.g.:
        %NOZZLES.specs{5}.D_or        = 5;
        %NOZZLES.specs{5}.V_ref       = 2;
        %NOZZLES.specs{5}.Delta_p_ref = 3*10^5;
        %NOZZLES.specs{5}.gamma_nz    = 30;

        % extraction from database only for not- or user- defined nozzles
        for i_nz = 1:length(NOZZLES.theta_nz)
            if ~( strcmp(NOZZLES.name_nz{i_nz},'un') || strcmp(NOZZLES.name_nz{i_nz},'ud'))
                [NOZZLES.specs{i_nz},fOK]=SX_DatabaseLoad('Nozzle',NOZZLES.name_nz{i_nz},fOK);
            end
        end, clear i_nz

        %% STRESS %%
        % Creation of STRESS structure: mechanical properties of CAST IRON
        STRESS.sigmaR_traz  = 250;  % [MPa]  Ultimate tensile stress
        STRESS.sigmaR_comp  = 850;  % [MPa]  Ultimate compressive stress
        STRESS.ratio_FA     = 0.25; % [adim] Fatigue stress ratio
        STRESS.b2           = 1;    % [adim] Fatigue size factor
        STRESS.b3           = 0.96; % [adim] Fatigue surface finish factor
        
        %% STRESS ON SHAFT %%
        % ================================================================
        % ================================================================
        % Creation of STRESSsh structure: mechanical properties of spheroidal cast iron (white cast iron) 
        STRESSsh.Rmt_s       = 450;    % [MPa] Ultimate tensile strength 
        STRESSsh.Rmc_s       = 700;    % [MPa] Ultimate compression strength 
        STRESSsh.Rp02_s      = 310;    % [MPa] Yield strength 
        STRESSsh.sigma_fas   = 210;    % [MPa] Fatigue strength 
        STRESSsh.b2_s        = 0.9;    % [adim] Dimension factor (recommended value: 0.9 if dhub > 0.015 [m] otherwise 1)            
        STRESSsh.b3_s        = 0.975;  % [adim] Surface roughness factor                    
        % ================================================================
        % ================================================================

        %% RUN SIMULATION %%
        % in case of new simulation, preprocessing and main are called
        if fOK
            SVECpreProc(IO, FLAG, NUMERIC, PROCESS, GEOMETRY, GAS, LIQUID, LEAK, VANE, NOZZLES, STRESS, STRESSsh, SHAFT, BUSHING);
        end

        % ================================================================
        % ================================================================
        % Addition of new structures as input for SVECpreProc
        % ================================================================
        % ================================================================
        
        if fOK
            [fOK] = SVECmodelMain(IO);
        end
    end
    
    % in case of saved simulation, only post-processing in called
    if fOK
        SVECpostProc(IO);
    end

    if fOK
        SX_Logfile('v',{'CALL_SVEC_M135A completed. Elapsed time: %.2f s ',toc(IO.CALLtime)})
    else
        fprintf('\nCALL_SVEC: an error occoured. SVEC could not complete the simulation. Read the previous warning messages or the logfile to try to fix the problem\n\n');
        SX_Logfile('a',{'Simulation aborted. Check the log file to try to fix the problem'})
    end
end
