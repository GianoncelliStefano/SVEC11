function SVECpreProc (IO, FLAG, NUMERIC, PROCESS, GEOMETRY, GAS, LIQUID, LEAK, VANE, NOZZLES, STRESS, STRESSsh, SHAFT, BUSHING,INTAKE)
% This function take as input all the structure from CALL_SVECmodelMain and
% perform some preprocessing operations:
% - conversion of string variables to numeric (rot_dir)
% - conversion to SI
% - grid discretization points correction (only if stress analysis is activated)
% - NOZZLES struct cleaning from deactivated nozzles
% - NOZZLES struct sorting
% Then, it saves a .mat files with the modified simulation parameters
%
% INPUT
% IO        : Structure of basic input output parameters
% FLAG      : Structure of flag parameters
% NUMERIC   : Structure of numeric parameters
% PROCESS   : Structure of process parameters
% GEOMETRY  : Structure of geometric parameters
% GAS       : Structure of gas parameters
% LIQUID    : Structure of oil parameters
% LEAK      : Structure of leakages parameters
% VANE      : Structure of vane parameters
% NOZZLES   : Structure of nozzles parameters
% STRESS    : Structure of stress parameters
% INTAKE    : Structure of intake model parameters
% =======================================================================
% =======================================================================
% STRESSsh  : Structure of shaft stress parameters
% SHAFT     : Structure of shaft parameters
% BUSHING   : Structure of bushing parameters
% =======================================================================
% =======================================================================
% 
%
% OUTPUT
% .mat file :this function generates a .mat file of the current simulation, containing
%            all the simulation parameters

    %% PREAMBLE %%
    PREtime = tic;
    SX_Logfile('v',{'PREprocessing started'})

    %% CONVERSION OF STRING VARIABLE TO NUMERIC %%
    % Conversion of the variable representing the rotational direction to:
    % - +1 clockwise;
    % - -1 counterclockwise;
    switch PROCESS.rot_dir
        case 'clockwise'
            PROCESS.rot_dir = 1;
        case 'counterclockwise'
            PROCESS.rot_dir = -1;
    end
    
    %% CONVERSION TO SI %%
    % Numeric struct
    NUMERIC.stp    = NUMERIC.stp *1e-3;  % [mm -> m] vane discretization step for stress analysis

    % Process struct
    PROCESS.p_suc = PROCESS.p_suc *1e+5;     % [bar -> Pa]  suction pressure
    PROCESS.p_del = PROCESS.p_del *1e+5;     % [bar -> Pa]  delivery pressure
    PROCESS.T_suc = PROCESS.T_suc + 273.15;  % [°C -> K]    gas inlet temperature
    PROCESS.T_0   = PROCESS.T_0 + 273.15;    % [°C -> K]    reference temperature for energy balance
    PROCESS.p_0   = PROCESS.p_0 *1e+5;       % [bar -> Pa]  reference pressure     
    PROCESS.zeta  = deg2rad(PROCESS.zeta);   % [deg -> rad] angles of x axis with respect to X

    % geometry struct
    c            = GEOMETRY.c;        % 0 circular stator;       1 elliptical stator;
    if c == 1
        GEOMETRY.D            = GEOMETRY.D *1e-3;               % [mm -> m] stator diameter
        GEOMETRY.d            = GEOMETRY.d *1e-3;               % [mm -> m] rotor diameter
        GEOMETRY.L            = GEOMETRY.L *1e-3;               % [mm -> m] rotor length
        GEOMETRY.l            = GEOMETRY.l *1e-3;               % [mm -> m] vane height
        GEOMETRY.s            = GEOMETRY.s *1e-3;               % [mm -> m] vane thickness
        GEOMETRY.s_1D         = GEOMETRY.s_1D *1e-3;            % [mm -> m] vane thickness for 1D geometry
        GEOMETRY.d_hub        = GEOMETRY.d_hub *1e-3;           % [mm -> m] bushing diameter
        GEOMETRY.th_SucOpen   = deg2rad(GEOMETRY.th_SucOpen);   % [deg -> rad] suction opening angle
        GEOMETRY.th_SucClose  = deg2rad(GEOMETRY.th_SucClose);  % [deg -> rad] suction close angle
        GEOMETRY.th_DisOpen   = deg2rad(GEOMETRY.th_DisOpen);   % [deg -> rad] delivery opening angle
        GEOMETRY.th_DisClose  = deg2rad(GEOMETRY.th_DisClose);  % [deg -> rad] delivery closing angle
        GEOMETRY.th_tilt      = deg2rad(GEOMETRY.th_tilt);      % [deg -> rad] vane tilt angle
        GEOMETRY.b            = GEOMETRY.b *1e-3;               % [mm -> m] offset between center of tip circonference and vane axis
        GEOMETRY.r_tip        = GEOMETRY.r_tip *1e-3;           % [mm -> m] vane tip radius
        GEOMETRY.TgAngle      = deg2rad(GEOMETRY.TgAngle);      % [°]  sealing arc between stator and rotor
        GEOMETRY.RSclr        = GEOMETRY.RSclr*1e-6;            % [micron->m] clearance between rotor and stator
        GEOMETRY.VSclr        = GEOMETRY.VSclr*1e-6;            % [micron->m] clearance between vane and end-plates
        GEOMETRY.PEclr        = GEOMETRY.PEclr*1e-6;            % [micron->m] clearance between rotor and end-plates
        GEOMETRY.INport_Amin  =  GEOMETRY.INport_Amin*1e-6;     % [mm^2-->m^2] Inlet port minimum passage area 
        GEOMETRY.INport_Amax  =  GEOMETRY.INport_Amax*1e-6;     % [mm^2-->m^2] Inlet port maximum passage area 
        GEOMETRY.OUTport_Amin =  GEOMETRY.OUTport_Amin*1e-6;    % [mm^2-->m^2] Outlet port minimum passage area
        GEOMETRY.OUTport_Amax =  GEOMETRY.OUTport_Amax*1e-6;    % [mm^2-->m^2] Outlet port maximum passage area 
    elseif c == 2
        GEOMETRY.d            = GEOMETRY.d *1e-3;               % [mm -> m] rotor diameter
        GEOMETRY.L            = GEOMETRY.L *1e-3;               % [mm -> m] rotor length
        GEOMETRY.l            = GEOMETRY.l *1e-3;               % [mm -> m] vane height
        GEOMETRY.s            = GEOMETRY.s *1e-3;               % [mm -> m] vane thickness
        GEOMETRY.d_hub        = GEOMETRY.d_hub *1e-3;           % [mm -> m] bushing diameter
        GEOMETRY.th_SucOpen   = deg2rad(GEOMETRY.th_SucOpen);   % [deg -> rad] suction opening angle
        GEOMETRY.th_SucClose  = deg2rad(GEOMETRY.th_SucClose);  % [deg -> rad] suction close angle
        GEOMETRY.th_DisOpen   = deg2rad(GEOMETRY.th_DisOpen);   % [deg -> rad] delivery opening angle
        GEOMETRY.th_DisClose  = deg2rad(GEOMETRY.th_DisClose);  % [deg -> rad] delivery closing angle
        GEOMETRY.TgAngle      = deg2rad(GEOMETRY.TgAngle);      % [°]  sealing arc between stator and rotor
        GEOMETRY.RSclr        = GEOMETRY.RSclr*1e-6;            % [micron->m] clearance between rotor and stator
        GEOMETRY.VSclr        = GEOMETRY.VSclr*1e-6;            % [micron->m] clearance between vane and end-plates
        GEOMETRY.PEclr        = GEOMETRY.PEclr*1e-6;            % [micron->m] clearance between rotor and end-plates
        GEOMETRY.INport_Amin  = GEOMETRY.INport_Amin*1e-6;      % [mm^2-->m^2] inlet port minimum passage area 
        GEOMETRY.INport_Amax  = GEOMETRY.INport_Amax*1e-6;      % [mm^2-->m^2] inlet port maximum passage area 
        GEOMETRY.OUTport_Amin = GEOMETRY.OUTport_Amin*1e-6;     % [mm^2-->m^2] outlet port minimum passage area
        GEOMETRY.OUTport_Amax = GEOMETRY.OUTport_Amax*1e-6;     % [mm^2-->m^2] outlet port maximum passage area 
    
    end, clear c
    
    if FLAG.fSDP
    % Intake struct
       INTAKE.lenght     = INTAKE.lenght*1e-3;                  % pipe length [m] 
       INTAKE.D_up       = INTAKE.D_up*1e-3;                    % diameter at pipe's start (upstream side) [m]
       INTAKE.D_do       = INTAKE.D_do*1e-3;                    % diameter at pipe's end (downstream side) [m]
       INTAKE.roughness  = INTAKE.roughness*1e-6;               % roughness [m]    
       INTAKE.cpitch     = INTAKE.cpitch*1e-3;                  % corrugated pitch [m]
       INTAKE.ct         = INTAKE.ct*1e-3;                      % corrugation height [m]
    end
    
    % Gas struct
    GAS.MM_g = GAS.MM_g*1e-3;                                   % [kg/kmol -> kg/mol] gas molar mass

    % Liquid struct (oil)
    LIQUID.T_ref    = LIQUID.T_ref +273.15;                     % [°C -> K]    oil reference temperature
    
    % =======================================================================
    % =======================================================================
    % Bushing struct
    BUSHING.d1 = BUSHING.d1*1e-3;     % [mm -> m] Rotor end-bushing distance 
    BUSHING.bB = BUSHING.bB*1e-3;     % [mm -> m] Bushing length 
    % =======================================================================
    % =======================================================================

    % Nozzles struct
    NOZZLES.p_inj    = NOZZLES.p_inj .*1e+5;          % [bar -> Pa] injection pressure
    NOZZLES.Tl_in    = NOZZLES.Tl_in + 273.15;        % [°C -> K]   oil inlet temperature temperature
    NOZZLES.theta_nz = deg2rad(NOZZLES.theta_nz);     % [deg -> rad] nozzle angular position

    % mechanical struct
    STRESS.sigmaR_traz = STRESS.sigmaR_traz *1e+6;    % [Mpa -> Pa] ultimate tensile stress
    STRESS.sigmaR_comp = STRESS.sigmaR_comp *1e+6;    % [Mpa -> Pa] ultimate compressive stress
    
    % =======================================================================
    % =======================================================================
    STRESSsh.Rmt_s      = STRESSsh.Rmt_s*1e+6;        % [MPa -> Pa] Ultimate tensile strength 
    STRESSsh.Rmc_s      = STRESSsh.Rmc_s*1e+6;        % [MPa -> Pa] Ultimate compression strength 
    STRESSsh.Rp02_s     = STRESSsh.Rp02_s*1e+6;       % [MPa -> Pa] Yield strength 
    STRESSsh.sigma_fas  = STRESSsh.sigma_fas*1e+6;    % [MPa -> Pa] Fatigue strength 
    % =======================================================================
    % =======================================================================
    %% GRID DISCRETIZATION POINTS CORRECTION %%
    % =====================================================================
    % =====================================================================
    if FLAG.fSTR == 1 && NUMERIC.Npt_i ~= 99000        
        warning ('S1_InputCheck:NUMERIC','Number of grid discretization point not compatible with Stress Analysis resampling method. SVEC will use Npt_i = 99000.');
        SX_Logfile ('w',{lastwarn});
        NUMERIC.Npt_i = 99000;
    end
    % =====================================================================
    % =====================================================================

    %% NOZZLE CLEANING %%
    % NOZZLES stucture is cleaned from deactiveted nozzles
    cnd                   = find(NOZZLES.f_nz==0);
    if ~isempty (cnd)
        warning('SVECpreprocessing:theta','Nozzles number %s have been deactivated, according to user input', mat2str(cnd'));
        SX_Logfile('v',{lastwarn});
        NOZZLES.theta_nz(cnd) = [];
        NOZZLES.num_nz(cnd)   = [];
        NOZZLES.p_inj(cnd)    = [];
        NOZZLES.Tl_in(cnd)    = [];
        NOZZLES.num_cls(cnd)  = [];
        NOZZLES.alpha_nz(cnd) = [];
        NOZZLES.type_nz(cnd)  = [];
        NOZZLES.name_nz(cnd)  = [];
        NOZZLES.specs(cnd)    = [];
        NOZZLES.f_nz (cnd)    = [];
    end

    %% NOOZLES SORTING %%
    % Nozzles are sorted in ascending order
    if ~issorted (NOZZLES.theta_nz, 'ascend')
        warning('SVECpreprocessing:theta_nz','NOZZLES.theta_nz must be in sorted ascending order. SVEC will sort the nozzles for you');
        SX_Logfile('v',{lastwarn})
        [~,idx]          = sort(NOZZLES.theta_nz, 'ascend');
        NOZZLES.f_nz     = NOZZLES.f_nz(idx);
        NOZZLES.theta_nz = NOZZLES.theta_nz(idx);
        NOZZLES.num_nz   = NOZZLES.num_nz(idx);
        NOZZLES.p_inj    = NOZZLES.p_inj(idx);
        NOZZLES.Tl_in    = NOZZLES.Tl_in(idx);
        NOZZLES.num_cls  = NOZZLES.num_cls(idx);
        NOZZLES.alpha_nz = NOZZLES.alpha_nz(idx);
        NOZZLES.type_nz  = NOZZLES.type_nz(idx);
        NOZZLES.name_nz  = NOZZLES.name_nz(idx);
        NOZZLES.specs    = NOZZLES.specs(idx);
    end, clear idx cnd

    %% SAVE SIMULATION FILE %%
    % save simulation.mat file except for IO struct  
    save (IO.fullsavename, 'FLAG', 'NUMERIC', 'PROCESS', 'GEOMETRY', 'GAS', 'LIQUID', 'LEAK', 'VANE', 'NOZZLES', 'STRESS', 'STRESSsh', 'SHAFT', 'BUSHING','INTAKE');
    SX_Logfile('v',{'PREprocessing completed. Elapsed time: %.2f s ',toc(PREtime)});
    
    % ================================================================
    % ================================================================
    % Addition of new structures as input for save
    % ================================================================
    % ================================================================
        
end
