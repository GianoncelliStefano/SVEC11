function [fOK] = S1_InputCheck(IO,FLAG, NUMERIC, PROCESS, GEOMETRY, GAS, LEAK, VANE, NOZZLES)
% This function take as input all the structure from CALL_SVECmodelMain and
% check if the input from user are correct and consistent. If an error or uncoherence 
% is spotted, a warning message is shown, the SVEC will exit the simulation
% INPUT
% IO        : Structure of basic input output parameters
% FLAG      : Structure of flag parameters
% NUMERIC   : Structure of numeric parameters
% GEOMETRY  : Structure of geometric parameters
% GAS       : Structure of gas parameters
% LEAK      : Structure of leakage parameters
% VANE      : Structure of vane parameters
% NOZZLES   : Structure of nozzles parameters
%
%OUTPUT
% fOK      : flag for input structure
%             0: structures not correclty implemented, exit the simulation
%             1: structures correclty implemented

    %% DEFINITION %%
    fOK = 1;
    
   %% IO INPUT CHECK %%
    f1  = IO.fMODE ~= 0 && IO.fMODE ~= 1;
    f2  = IO.fDPLT ~= 0 && IO.fDPLT ~= 1;
    f3  = IO.fDRPT ~= 0 && IO.fDRPT ~= 1;
    f4  = IO.fSRPT ~= 0 && IO.fSRPT ~= 1; 
    f5  = IO.fSPLT ~= 0 && IO.fSPLT ~= 1;
    
    if f1
        warning('S1_InputCheck:IO','Unknown value of IO.fMODE. IO variables can either be 0 or 1');
        SX_Logfile ('e',{lastwarn});
        fOK = 0;
    end
    
    if f2
        warning('S1_InputCheck:IO','Unknown value of IO.fDPLT. IO variables can either be 0 or 1');
        SX_Logfile ('e',{lastwarn});
        fOK = 0;
    end
    if f3
        warning('S1_InputCheck:IO','Unknown value of IO.fDRPT. IO variables can either be 0 or 1');
        SX_Logfile ('e',{lastwarn});
        fOK = 0;
    end
    if f4
        warning('S1_InputCheck:IO','Unknown value of IO.fSPLT. IO variables can either be 0 or 1');
        SX_Logfile ('e',{lastwarn});
        fOK = 0;
    end
    if f5
        warning('S1_InputCheck:IO','Unknown value of IO.fSRPT. IO variables can either be 0 or 1');
        SX_Logfile ('e',{lastwarn});
        fOK = 0;
    end
    clear f1 f2 f3 f4 f5
    
    %% FLAG INPUT CHECK %%
    f1  = FLAG.fDBG ~= 0 && FLAG.fDBG ~= 1;
    f2  = FLAG.fSTR ~= 0 && FLAG.fSTR ~= 1;
    f3  = FLAG.fSDP ~= 0 && FLAG.fSDP ~= 1;
    if f1
        warning('S1_InputCheck:FLAG','Unknown value of FLAG.fDBG. FLAG variables can either be 0 or 1');
        SX_Logfile ('e',{lastwarn});
        fOK = 0;
    end
    if f2
        warning('S1_InputCheck:FLAG','Unknown value of FLAG.fSTR. FLAG variables can either be 0 or 1');
        SX_Logfile ('e',{lastwarn});
        fOK = 0;
    end
    if f3
        warning('S1_InputCheck:FLAG','Unknown value of FLAG.fSDP. FLAG variables can either be 0 or 1');
        SX_Logfile ('e',{lastwarn});
        fOK = 0;
    end
    clear f1 f2 f3
    
    %% NUMERIC INPUT CHECK %%
    % =====================================================================
    % =====================================================================
    if FLAG.fSTR == 1
       f1 = NUMERIC.Npt_i ~= 99000;
       
       if f1
           NUMERIC.Npt_i = 99000;
           warning ('S1_InputCheck:NUMERIC','Number of grid discretization point not compatible with Stress Analysis resampling method. Npt_i has been imposed equal to 99000');
           SX_Logfile ('w',{lastwarn});
       end
    end
    clear f1
    % =====================================================================
    % =====================================================================
    %% PROCESS INPUT CHECK %%
    if PROCESS.process == 1 && (PROCESS.p_del <= PROCESS.p_suc)
        warning ('S1_InputCheck:PROCESS','Compression - Delivery pressure should be higher than suction pressure');
        SX_Logfile ('w',{lastwarn});
        
    elseif PROCESS.process == 2 && (PROCESS.p_del >= PROCESS.p_suc)
        warning ('S1_InputCheck:PROCESS','Expansion - Delivery pressure should be lower than suction pressure');
        SX_Logfile ('w',{lastwarn});
        
    elseif PROCESS.process ~= 1 && PROCESS.process ~= 2
        warning ('S1_InputCheck:PROCESS','Unknown value of PROCESS.process');
        SX_Logfile ('e',{lastwarn});
        fOK = 0;
    end
    
    %% GEOMETRY INPUT CHECK %%
    if GEOMETRY.c == 1                                              % circular geometry
        f1 = GEOMETRY.d > GEOMETRY.D;                               % rotor diameter must be less then stator diameter
        f2 = GEOMETRY.r_tip < 0;                                    % tip radius must be positive
        f3 = GEOMETRY.r_tip < 0.5*GEOMETRY.s + abs(GEOMETRY.b);     % limit condition for tip radius
        f4 = GEOMETRY.th_tilt*GEOMETRY.b < 0;                       % cuspidal point slippage may happen - tilt angle and tip offset better have same sign
        f5 = GEOMETRY.l > GEOMETRY.d;                               % vane height must be less than rotor diameter
        
        if f1
            warning('S1_InputCheck:GEOMETRY','Rotor diameter is greater then stator diameter');
            SX_Logfile ('e',{lastwarn});
            fOK = 0;
        end
        if f2
            warning('S1_InputCheck:GEOMETRY','r_tip is negative');
            SX_Logfile ('e',{lastwarn});
            fOK = 0;
        end
        if f3
            warning('S1_InputCheck:GEOMETRY','Unfeasable vane tip geoemtry: please modify b and r_tip');
            SX_Logfile ('e',{lastwarn});
            fOK = 0;
        end
        if f4
            warning('S1_InputCheck:GEOMETRY','Theta_tilt and b have opposite signs: cuspidale point slippage may happen');
            SX_Logfile ('w',{lastwarn});
        end
        if f5
            warning('S1_InputCheck:GEOMETRY','Vane length greater then rotor diameter');
            SX_Logfile ('e',{lastwarn});
            fOK = 0;
        end

    elseif GEOMETRY.c == 2                    % elliptical
        f6 = GEOMETRY.e <0 || GEOMETRY.e >1;  % eccentricity must be between 0 and 1
        f7 = GEOMETRY.l > GEOMETRY.d;         % vane height must be less than rotor diameter

        if f6
            warning ('S1_InputCheck:GEOMETRY','Eccentricity ''e'' must be between 0 and 1')
            SX_Logfile ('e',{lastwarn});
            fOK = 0;
        end
        if f7
            warning('S1_InputCheck:GEOMETRY','Vane length greater then rotor diameter');
            SX_Logfile ('e',{lastwarn});
            fOK = 0;
        end
    else
        warning ('S1_InputCheck:LOGIC','Unknown value of GEOMETRY.c. Please select 0 or 1');
        SX_Logfile ('e',{lastwarn});
        fOK = 0;
    end

    % check if angles position and order is ok
    f8  = GEOMETRY.th_SucOpen  > (2*pi/GEOMETRY.c) || GEOMETRY.th_SucOpen  <0;
    f9  = GEOMETRY.th_SucClose > (2*pi/GEOMETRY.c) || GEOMETRY.th_SucClose <0;
    f10 = GEOMETRY.th_DisOpen  > (2*pi/GEOMETRY.c) || GEOMETRY.th_DisOpen  <0;
    f11 = GEOMETRY.th_DisClose > (2*pi/GEOMETRY.c) || GEOMETRY.th_DisClose <0;

    if (f8 || f9 || f10 || f11) == 1
        warning ('S1_InputCheck:GEOMETRY','Ports angles are not included into 0:%i ° interval',360/GEOMETRY.c)
        SX_Logfile ('e',{lastwarn});
        fOK = 0;
    end
    
    % check design process
    f12 = GEOMETRY.StdProcess ~= 1 && GEOMETRY.StdProcess ~= 2;  % Design process can either be compression or expansion
    if f12
        warning('S1_InputCheck:GEOMETRY','Design process can either be 1 (compression) or 2 (expansion)');
        SX_Logfile ('e',{lastwarn});
        fOK = 0;
    end
    
    % check volume maximization (angle modifier)
    f13 = (GEOMETRY.StdProcess ~= PROCESS.process) && ( isnan(GEOMETRY.th_SucClose) || isnan(GEOMETRY.th_DisOpen) );
    f14 = isnan(GEOMETRY.th_SucClose) && PROCESS.process == 2;
    f15 = isnan(GEOMETRY.th_DisOpen) && PROCESS.process == 1;
    
    if f13
        warning('S1_InputCheck:GEOMETRY','Volume maximization done by ports change is available only if design process and performed process match');
        SX_Logfile ('e',{lastwarn});
        fOK = 0;
    end
    
    if f14
        warning('S1_InputCheck:GEOMETRY','Volume maximization done by suction close change is available only with COMPRESSION process ');
        SX_Logfile ('e',{lastwarn});
        fOK = 0;
    end
    
    if f15
        warning('S1_InputCheck:GEOMETRY','Volume maximization done by delivery close change is available only with EXPANSION process ');
        SX_Logfile ('e',{lastwarn});
        fOK = 0;
    end
    
    % check suction and discharge port areas 
    f16 = GEOMETRY.INport_Amax  < GEOMETRY.INport_Amin;
    f17 = GEOMETRY.OUTport_Amax < GEOMETRY.OUTport_Amin;
    
    if f16
        warning('S1_InputCheck:GEOMETRY','Inlet port maximum area is intended to be larger than the minimum area');
        SX_Logfile ('e',{lastwarn});
        fOK = 0;
    end
    
    if f17
        warning('S1_InputCheck:GEOMETRY','Outlet port maximum area is intended to be larger than the minimum area');
        SX_Logfile ('e',{lastwarn});
        fOK = 0;
    end
    
    clear f1 f2 f3 f4 f5 f6 f7 f8 f9 f10 f11 f12 f13 f14 f15 f16 f17
    
    %% GAS INPUT CHECK %%
    f1 = abs(sum(GAS.molcomp_g)-1)>NUMERIC.toll_d;                               % sum of molar fraction must be one
    f2 = ~strcmp(GAS.model_g, 'simple') && ~strcmp(GAS.model_g,'full');          % model selection is limited
    f3 = ~strcmp(GAS.propSource, 'direct') && ~strcmp(GAS.propSource,'REFPROP') && ~strcmp(GAS.propSource,'gas'); % source of thermodynamic data is limited
    f4 = ~strcmp(GAS.model_full, 'fullTV') && ~strcmp(GAS.model_full, 'fullTP'); % model selection is limited
    
    if f1
        warning ('S1_InputCheck:GAS','gas molar composition GAS.z is not correct')
        SX_Logfile ('e',{lastwarn});
        fOK = 0;
    end
    if f2
        warning ('S1_InputCheck:GAS','GAS.model not found. Plese select ''simple'' or ''full'' ')
        SX_Logfile ('w',{lastwarn});
        fOK = 0;
    end
    if f3
        warning ('S1_InputCheck:GAS','GAS.propSource not found. Plese select ''direct'' or ''REFPROP'' ')
        SX_Logfile ('w',{lastwarn});
        fOK = 0;
    end
    if f4
        warning ('S1_InputCheck:GAS','GAS.model_full not found. Plese select ''fullTV'' or ''fullTP'' ')
        SX_Logfile ('w',{lastwarn});
        fOK = 0;
    end
    clear f1 f2 f3 f4
    
    %% LEAK INPUT CHECK %%
    % Rotor-stator
    f1 = ~strcmp(LEAK.MUmodel_RS, 'Dukler') && ~strcmp(LEAK.MUmodel_RS,'Maxwell') && ~strcmp(LEAK.MUmodel_RS,'Awad') && ~strcmp(LEAK.MUmodel_RS,'MacAdams') && ~strcmp(LEAK.MUmodel_RS,'Cicchitti') && ~strcmp(LEAK.MUmodel_RS,'PureOil');                                 % model selection is limited
    f2 = ~strcmp(LEAK.RSmodel, 'Ishii') && ~strcmp(LEAK.RSmodel,'Ferreira') && ~strcmp(LEAK.RSmodel,'Poiselle-Couette') && ~strcmp(LEAK.RSmodel,'Gasche') && ~strcmp(LEAK.RSmodel,'Yanagisawa') && ~strcmp(LEAK.RSmodel,'Badr');  % model selection is limited
    f3 = strcmp(LEAK.MUmodel_RS, 'PureOil') && sum(NOZZLES.f_nz==0);
    % Vane-side
    f4 = ~strcmp(LEAK.MUmodel_VS, 'Dukler') && ~strcmp(LEAK.MUmodel_VS,'Maxwell') && ~strcmp(LEAK.MUmodel_VS,'Awad') && ~strcmp(LEAK.MUmodel_VS,'MacAdams') && ~strcmp(LEAK.MUmodel_VS,'Cicchitti') && ~strcmp(LEAK.MUmodel_VS,'PureOil');                                 % model selection is limited
    f5 = ~strcmp(LEAK.VSmodel, 'Ishii') && ~strcmp(LEAK.VSmodel,'Yuan') && ~strcmp(LEAK.VSmodel,'Poiselle-Couette') && ~strcmp(LEAK.VSmodel,'Yanagisawa') && ~strcmp(LEAK.VSmodel,'Badr') && ~strcmp(LEAK.VSmodel,'Suefuji');  % model selection is limited
    f6 = strcmp(LEAK.MUmodel_VS, 'PureOil') && sum(NOZZLES.f_nz==0);
    % Rotor-end-plate
    f7 = ~strcmp(LEAK.MUmodel_PE, 'Dukler') && ~strcmp(LEAK.MUmodel_PE,'Maxwell') && ~strcmp(LEAK.MUmodel_PE,'Awad') && ~strcmp(LEAK.MUmodel_PE,'MacAdams') && ~strcmp(LEAK.MUmodel_PE,'Cicchitti') && ~strcmp(LEAK.MUmodel_PE,'PureOil');                                 % model selection is limited
    f8 = ~strcmp(LEAK.PEmodel, 'Badr') && ~strcmp(LEAK.PEmodel,'Poiselle-Couette') && ~strcmp(LEAK.PEmodel, 'Ishii');  % model selection is limited
    f9 = strcmp(LEAK.MUmodel_PE, 'PureOil') && sum(NOZZLES.f_nz==0);
    
    if f1
        warning ('S1_InputCheck:LEAK','LEAK.MUmodel_RS: Select a proper mixture viscosity model')
        SX_Logfile ('e',{lastwarn});
        fOK = 0;
    end
    
    if f2
        warning ('S1_InputCheck:LEAK','LEAK.RSmodel: Select a proper rotor-stator leakage model')
        SX_Logfile ('e',{lastwarn});
        fOK = 0;
    end
    
    if f3
        warning ('S1_InputCheck:LEAK','LEAK.MUmodel_RS: Cannot use ''PureOil model'' if nozzle are all disabled')
        SX_Logfile ('e',{lastwarn});
        fOK = 0;
    end
    
    if f4
        warning ('S1_InputCheck:LEAK','LEAK.MUmodel_RS: Select a proper mixture viscosity model')
        SX_Logfile ('e',{lastwarn});
        fOK = 0;
    end
    
    if f5
        warning ('S1_InputCheck:LEAK','LEAK.VSmodel: Select a proper vane-side leakage model')
        SX_Logfile ('e',{lastwarn});
        fOK = 0;
    end
    
    if f6
        warning ('S1_InputCheck:LEAK','LEAK.MUmodel_VS: Cannot use ''PureOil model'' if nozzle are all disabled')
        SX_Logfile ('e',{lastwarn});
        fOK = 0;
    end
    
    if f7
        warning ('S1_InputCheck:LEAK','LEAK.MUmodel_PE: Select a proper mixture viscosity model')
        SX_Logfile ('e',{lastwarn});
        fOK = 0;
    end
    
    if f8
        warning ('S1_InputCheck:LEAK','LEAK.PEmodel: Select a proper rotor-end-plate leakage model')
        SX_Logfile ('e',{lastwarn});
        fOK = 0;
    end
    
    if f9
        warning ('S1_InputCheck:LEAK','LEAK.MUmodel_PE: Cannot use ''PureOil model'' if nozzle are all disabled')
        SX_Logfile ('e',{lastwarn});
        fOK = 0;
    end
    
    % Two-phase viscosity models and leakages compatibility: Problems might occur when:
    % - using Poiselle-Couette leakage model coupled with Dukler viscosity model (almost pure gas)
    % - using Poiselle-Couette leakage model in case of oil free machine (whatever the two-phase viscosity model)
%     f10 = strcmp(LEAK.MUmodel_VS, 'Dukler') && strcmp(LEAK.VSmodel,'Poiselle-Couette');
%     f11 = strcmp(LEAK.VSmodel,'Poiselle-Couette') && sum(NOZZLES.f_nz) == 0;
%     
%     if f10
%         warning ('S1_InputCheck:LEAK','LEAK.MUmodel_VS & LEAK.VSmodel: Cannot use ''Dukler'' model coupled with ''Poiselle-Couette'' model. Bulk viscosity too low')
%         SX_Logfile ('e',{lastwarn});
%         fOK = 0;
%     end
%     
%     if f11
%         warning ('S1_InputCheck:LEAK','LEAK.VSmodel: Cannot use ''Poiselle-Couette'' model if nozzle are all disabled. Bulk viscosity too low')
%         SX_Logfile ('e',{lastwarn});
%         fOK = 0;
%     end
%     clear f1 f2 f3 f4 f5 f6 f7 f8 f9 f10 f11
    
    %% VANE INPUT CHECK %%
    f1 = VANE.offsetG >100 || VANE.offsetG <0;                % offset of mass center is given as percentage
    f2 =  ~strcmp(VANE.vane_type, 'Standard') && FLAG.fSTR == 1;   % stress analysis check
    
    if f1
        warning ('S1_InputCheck:VANE','VANE.offsetG: Mass center of vane falls outside the vane')
        SX_Logfile ('e',{lastwarn});
        fOK = 0;
    end
    if f2
        warning('S1_InputCheck:VANE','Stress analysis is valid only for ''Standard'' vanes');
        SX_Logfile ('w',{lastwarn});
    end
    clear f1 f2
    
    %% NOZZLES CHECK %%   
    % Check on dimension of nozzles structure field
     if length(NOZZLES.f_nz) ~= length(NOZZLES.theta_nz)                
        warning('S1_InputCheck:NOZZLES_f_nz','NOZZLES.f_nz and NOZZLES.theta_nz must have same dimensions');
        SX_Logfile ('e',{lastwarn});
        fOK = 0;
    end
    if length(NOZZLES.num_nz) ~= length(NOZZLES.theta_nz)                
        warning('S1_InputCheck:NOZZLES_num_nz','NOZZLES.num_nz and NOZZLES.theta_nz must have same dimensions');
        SX_Logfile ('e',{lastwarn});
        fOK = 0;
    end
    if length(NOZZLES.p_inj) ~= length(NOZZLES.theta_nz)
       warning('S1_InputCheck:NOZZLES_p_inj','NOZZLES.p_inj and NOZZLES.theta_nz must have same dimensions');
       SX_Logfile ('e',{lastwarn});
       fOK = 0;
    end
    if length(NOZZLES.Tl_in) ~= length(NOZZLES.theta_nz)
        warning('S1_InputCheck:NOZZLES_Tl_in','NOZZLES.Tl_in and NOZZLES.theta_nz must have same dimensions');
        SX_Logfile ('e',{lastwarn});
        fOK = 0;
    end
    if length(NOZZLES.num_cls) ~= length(NOZZLES.theta_nz)
        warning('S1_InputCheck:NOZZLES_numcls','NOZZLES.numcls and NOZZLES.theta_nz must have same dimensions');
        SX_Logfile ('e',{lastwarn});
        fOK = 0;
    end
    if length(NOZZLES.alpha_nz) ~= length(NOZZLES.theta_nz)
        warning('S1_InputCheck:NOZZLES_alpha_nz','NOZZLES.alpha_nz and NOZZLES.theta_nz must have same dimensions');
        SX_Logfile ('e',{lastwarn});
        fOK = 0;
    end
    if length(NOZZLES.type_nz) ~= length(NOZZLES.theta_nz)
        warning('S1_InputCheck:NOZZLES_type_nz','NOZZLES.type_nz and NOZZLES.theta_nz must have same dimensions');
        SX_Logfile ('e',{lastwarn});
        fOK = 0;
    end
    if length(NOZZLES.name_nz) ~= length(NOZZLES.theta_nz)
        warning('S1_InputCheck:NOZZLES_name_nz','NOZZLES.name_nz and NOZZLES.theta_nz must have same dimensions');
        SX_Logfile ('e',{lastwarn});
        fOK = 0;
    end   
    if length(NOZZLES.specs) ~= length(NOZZLES.theta_nz)
        warning('S1_InputCheck:NOZZLES_specs','NOZZLES.specs and NOZZLES.theta_nz must have same dimensions');
        SX_Logfile ('e',{lastwarn});
        fOK = 0;
    end
    
   
    for i_nz=1:length(NOZZLES.theta_nz)
        if NOZZLES.f_nz(i_nz) == 1
            
            % check on nozzles droplet frequency
            if NOZZLES.alpha_nz(i_nz) >1 || NOZZLES.alpha_nz(i_nz) <0
                warning('S1_InputCheck:NOZZLES_alpha_nz','Active nozzle (%i): alpha_nz must be between 0 and 1' , i_nz);
                SX_Logfile ('e',{lastwarn});
                fOK = 0;
            end
            
            % check on nozzles types for all nozzles
            if ~( strcmp(NOZZLES.name_nz{i_nz},'un') || strcmp(NOZZLES.name_nz{i_nz},'ud') )
                if ~strncmp(NOZZLES.type_nz{i_nz},NOZZLES.specs{i_nz}.type,2)
                    warning('S1_InputCheck:NOZZLES_name_nz','Active nozzle (%i): Wrong type of nozzle selected. SVEC will consider the nozzles type in database as the right one',i_nz);
                    NOZZLES.type_nz{i_nz} = NOZZLES.specs{i_nz}.type;
                    SX_Logfile ('w',{lastwarn});
                end
            end
            
            if strncmp(NOZZLES.type_nz{i_nz},'po',2)
                if ~(strcmp(NOZZLES.type_nz{i_nz},'po') || strcmp(NOZZLES.type_nz{i_nz},'po2') || strcmp(NOZZLES.type_nz{i_nz},'po3') || strcmp(NOZZLES.type_nz{i_nz},'po4'))
                    warning('S1_InputCheck:NOZZLES_type_nz','Active nozzle (%i): unknown type of nozzle. ''po'' nozzles can be po, po2, po3 or po4' , i_nz);
                    SX_Logfile ('e',{lastwarn});
                    fOK = 0;
                end
            elseif strncmp(NOZZLES.type_nz{i_nz},'fc',2)
                if ~(strcmp(NOZZLES.type_nz{i_nz},'fc') || strcmp(NOZZLES.type_nz{i_nz},'fc2') || strcmp(NOZZLES.type_nz{i_nz},'fc3'))
                    warning('S1_InputCheck:NOZZLES_type_nz','Active nozzle (%i): unknown type of nozzle. ''fc'' nozzles can be fc, fc2 or fc3' , i_nz);
                    SX_Logfile ('e',{lastwarn});
                    fOK = 0;
                end
            elseif ~(strcmp(NOZZLES.type_nz{i_nz},'hc') || strcmp(NOZZLES.type_nz{i_nz},'un') || strcmp(NOZZLES.type_nz{i_nz},'kf'))
                warning('S1_InputCheck:NOZZLES_type_nz','Active nozzle (%i): unknown type of nozzle. Nozzle can have type ''po'',''fc'',''hc'', ''kf'' or ''un'' '  , i_nz);
                SX_Logfile ('e',{lastwarn});
                fOK = 0;
            end
            
            % Check on nozzles specs field for ud and un nozzles
            if strcmp(NOZZLES.name_nz{i_nz},'un') == 1
                if sum(isfield(NOZZLES.specs{i_nz},{'D_g' 'f_g' 'eps_m'}))== 3
                    if sum(NOZZLES.specs{i_nz}.f_g) ~= 1
                        warning('S1_InputCheck:NOZZLES_specs','Active nozzle (%i): sum of NOZZLES.specs.f_g must be equal to 1 ',i_nz);
                        SX_Logfile ('e',{lastwarn});
                        fOK = 0;
                    end
                    if length(NOZZLES.specs{i_nz}.D_g) ~= NOZZLES.num_cls(i_nz)
                        warning('S1_InputCheck:NOZZLES_specs','Active nozzle (%i): number of elements in NOZZLES.specs.D_g must be equal to NOZZLES.numcls',i_nz);
                        SX_Logfile ('e',{lastwarn});
                        fOK = 0;
                    end
                    if length(NOZZLES.specs{i_nz}.f_g) ~= NOZZLES.num_cls(i_nz)
                        warning('S1_InputCheck:NOZZLES_specs','Active nozzle (%i): number of elements in NOZZLES.specs.f_g must be equal to NOZZLES.numcls',i_nz);
                        SX_Logfile ('e',{lastwarn});
                        fOK = 0;
                    end
                else
                    warning('S1_InputCheck:NOZZLES_specs','Active nozzle (%i): missing field in NOZZLES.specs{%i}. For ''un'' nozzles, NOZZLES.specs must contain 3 fields: D_g, f_g eps_m',i_nz, i_nz);
                    SX_Logfile ('e',{lastwarn});
                    fOK = 0;
                end
            elseif strcmp(NOZZLES.name_nz{i_nz},'ud') == 1
                if strncmp(NOZZLES.type_nz{i_nz},'kf',2)
                    if sum(isfield(NOZZLES.specs{i_nz},{'D_g' 'f_g' 'V_inj'}))~= 3
                        warning('S1_InputCheck:NOZZLES_specs','Active nozzle (%i): missing field in NOZZLES.specs{%i}. For userdefined ''kf'' nozzles, NOZZLES.specs must contain 3 fields: D_g, f_g, V_inj',i_nz, i_nz);
                        SX_Logfile ('e',{lastwarn});
                        fOK = 0;
                    end
            elseif strncmp(NOZZLES.type_nz{i_nz},'po',2)
                    if sum(isfield(NOZZLES.specs{i_nz},{'D_or' 'L_or'}))~= 2
                        warning('S1_InputCheck:NOZZLES_specs','Active nozzle (%i): missing field in NOZZLES.specs{%i}. For userdefined ''po'' nozzles, NOZZLES.specs must contain 2 fields: D_or, L_or',i_nz, i_nz);
                        SX_Logfile ('e',{lastwarn});
                        fOK = 0;
                    end
                    if NOZZLES.num_cls ~= 1
                        warning('S1_InputCheck:NOZZLES_specs','Active nozzle (%i): For ''po'' nozzles, the number of class must be equal to 1',i_nz);
                        SX_Logfile ('e',{lastwarn});
                        fOK = 0;
                    end
                elseif strncmp(NOZZLES.type_nz{i_nz},'hc',2)
                    if sum(isfield(NOZZLES.specs{i_nz},{'D_or' 'V_ref' 'Delta_p_ref' 'gamma'}))~= 4
                        warning('S1_InputCheck:NOZZLES_specs','Active nozzle (%i): missing field in NOZZLES.specs{%i}. For userdefined ''hc'' nozzles, NOZZLES.specs must contain 4 fields: D_or, V_ref, Deltap_ref, gamma', i_nz, i_nz);
                        SX_Logfile ('e',{lastwarn});
                        fOK = 0;
                    end
                elseif strncmp(NOZZLES.type_nz{i_nz},'fc',2)
                    if sum(isfield(NOZZLES.specs{i_nz},{'D_or' 'V_ref' 'Delta_p_ref'}))~= 3
                        warning('S1_InputCheck:NOZZLES_specs','Active nozzle (%i): missing field in NOZZLES.specs{%i}. For userdefined ''fc'' nozzles, NOZZLES.specs must contain 3 fields: D_or, V_ref, Deltap_ref', i_nz, i_nz);
                        SX_Logfile ('e',{lastwarn});
                        fOK = 0;
                    end
                end
            end
        end
    end

end

