function [pos_SucOpen,pos_SucClose,pos_DisOpen,pos_DisClose,pos_End,th_in,th_out,th_flp,th_Vmax,fOK] = S2_Angles (c,n_van,process,StdProcess,theta,th_SucOpen,th_SucClose,th_DisOpen,th_DisClose,TgAngle,V_cell,pos_tang,sgn,psi,Gamma,NOZZLES)
% This function perform a series of changes on the ports angle defined by user:
% - correction due to expansion process
% - correction due to logic of suction close
% - correction due to tilted vanes
% - correction due to discretization
%
% INPUT
% c                  : geometry parameter (1 cilindrical / 2 elliptical)
% n_van [-]          : number of vanes
% process [-]        : logic parameter (1 compression / 2 expansion)
% StdProcess [-]     : design process performed by the machine (1 compression / 2 expansion)
% StdProcess [-]     : design process for current machine (1 compression / 2 expansion)
% theta [rad]        : array of discretized angular position of trailing vane
% th_SucOpen [rad]   : suction open angle
% th_SucClose [rad]  : suction close angle
% th_DisOpen [rad]   : delivery open angle
% th_DisClose [rad]  : delivery close angle
% TgAngle [rad]      : sealing arc between stator and rotor
% V_cell [m3]        : volume of a cell in range [-gamma:2pi]
% pos_tang [-]       : tangency position in array theta
% sgn                : vane tilting direction (positive/zero/negative for forard/radial/backward)          
% psi [rad]          : angle BOC
% Gamma [rad]        : angular extention of a cell (cavity between 2 vanes)
% NOZZLES            : nozzles input parameters

% OUTPUT
% pos_Start [-]      : discretized position of cell start in array theta (is where tangency ends)
% pos_SucOpen [-]    : discretized suction close angle position in array theta
% pos_SucClose [-]   : discretized suction close angle position in array theta
% pos_DisOpen [-]    : discretized delivery open angle position in array theta
% pos_DisClose [-]   : discretized delivery close angle position in array theta
% pos_End [-]        : discretized position of cell end in array theta (is where tangency starts)
% th_in [rad]        : updated suction close angle (takes into account if user wants to maximise suction volume & expansion ports flipping)
% th_out [rad]       : updated delivery open angle (takes into account expansion ports flipping)
% th_flp [rad]       : flipped angles
% th_Vmax [rad]      : port angle that maximize cell volume
% fOK [bool]         : OK flag (1 - no problem occourred, 0 - a problem have been spotted, SVEC will exit the simulation)
%
% NOTE: There can be only one of these two correction:
%          * Ports flipping correction
%          * Volume Maximization correction
%       In fact, if ports are flipped, the suction close/ delivery open
%       cannot be changed, since it is already known. In S1_InputCheck it
%       is verified that only one correction takes place

    %% PREAMBLE %%
    fOK      = 1;               % error flag  
    th_Start = TgAngle/2;       % cell start angle [rad]
    th_End   = 2*pi-TgAngle/2;  % cell end angle   [rad]
    
    %% PORTS FLIPPING CORRECTION %%
    % if the design process and the performed process do not match,
    % ports flipping is needed
    if process ~= StdProcess
        th_Start_flp    = 2*pi/c - th_End;
        th_SucOpen_flp  = 2*pi/c - th_DisClose;
        th_SucClose_flp = 2*pi/c - th_DisOpen;
        th_DisOpen_flp  = 2*pi/c - th_SucClose;
        th_DisClose_flp = 2*pi/c - th_SucOpen;
        th_End_flp      = 2*pi/c - th_Start;
        % ports angle must be renamed here in order to avoid conflict in calculations below
        th_Start    = th_Start_flp;
        th_SucOpen  = th_SucOpen_flp;
        th_SucClose = th_SucClose_flp;
        th_DisOpen  = th_DisOpen_flp;
        th_DisClose = th_DisClose_flp;
        th_End      = th_End_flp;
        th_flp = [th_SucOpen,th_SucClose,th_DisOpen,th_DisClose];
        warning ('S2_Angles:SucPort','Design process and performed process do not match. Ports angle flipping required')
        SX_Logfile ('v',{lastwarn});
    else
        th_flp = NaN;
    end
    
    %% VOLUME MAXIMIZATION CORRECTION %%
    % if volume maximization is needed, one angle needs a correction
    if isnan(th_SucClose) || isnan(th_DisOpen)
        [~,pos_maxV] = max(V_cell);   % angular position of max volume
        
        switch process
            case 1 % compression, change th_SucClose
                th_SucClose  = theta(pos_maxV);
                th_Vmax = th_SucClose;
                warning ('S2_Angles:SucPort','Suction close angle changed to maximise suction volume. New suction close angle: %.2f°', rad2deg(th_SucClose));
                SX_Logfile ('v',{lastwarn});
            case 2 % expansion, change th_DisOpen
                th_DisOpen  = 2*pi/c - theta(pos_maxV);
                th_Vmax = th_DisOpen;
                warning ('S2_Angles:DisOpen','Delivery open angle changed to maximise delivery volume. New delivery open angle: %.2f°', rad2deg(th_DisOpen));
                SX_Logfile ('v',{lastwarn});
        end
    else
        th_Vmax = NaN;
    end
    
    %% TILTED VANES CORRECTION %%
    % Computes the angular position of ports as seen by rotor, in case of tilted vanes
    if c == 1 
    thetaVane = theta(pos_tang:end);
    m = thetaVane + sgn*psi;
    
    % referece indexes for tilted vanes
    [~,idx_th_Start]  = min(abs(m-th_Start));
    [~,idx_SucOpen]   = min(abs(m-th_SucOpen));
    [~,idx_SucClose]  = min(abs(m-th_SucClose));
    [~,idx_DisOpen]   = min(abs(m-th_DisOpen));
    [~,idx_DisClose]  = min(abs(m-th_DisClose));
    [~,idx_Th_End]    = min(abs(m-th_End));
    
    % new ports angle are computed
    th_Start    = th_Start - sgn* psi(idx_th_Start);
    th_SucOpen  = th_SucOpen - sgn* psi(idx_SucOpen);
    th_SucClose = th_SucClose - sgn* psi(idx_SucClose);
    th_DisOpen  = th_DisOpen - sgn* psi(idx_DisOpen);
    th_DisClose = th_DisClose - sgn* psi(idx_DisClose);
    th_End      = th_End - sgn* psi(idx_Th_End);
    clear thetaPala m idx_SucOpen idx_SucClose idx_DisOpen idx_DisClose     
    end
    
    %% DISCRETIZATION CORRECTION %%
    % Due to discretization, ports position is chosen as the closest one in array theta
    %[~,pos_Start]     = min(abs(theta-th_Start));   % cell start discretization position
    [~,pos_SucOpen] = min(abs(theta-th_SucOpen));    % suction open discretized position (for the future...)
    [~,pos_SucClose]  = min(abs(theta-th_SucClose)); % suction close discretized position
    [~,pos_DisOpen]   = min(abs(theta-th_DisOpen));  % delivery open discretized position
    [~,pos_DisClose]  = min(abs(theta-th_DisClose)); % delivery close discretized position (for the future...)
    [~,pos_End]       = min(abs(theta-th_End));      % cell end discretization position
    
    %% CLOSED CELL ANGLES %%
    % this angle will be used for plotting
    th_in  = th_SucClose;
    th_out = th_DisOpen;
    
    %% CHECK %%
    % check that after the modification, angle position is still ok. It may
    % also happen that by choosing the suction close angle that maximise the
    % volume, this angle will be lower than the suction open angle
    
    f1 = ~issorted([th_SucOpen,th_SucClose,th_DisOpen,th_DisClose],'monotonic');
    f2 = (th_DisOpen - th_SucClose) < 2*pi/n_van;
    
    if f1 
        warning ('S2_Angles:PortsOrder','Ports order is not correct or ports angles are not correct')
        fOK = 0;
        SX_Logfile ('e',{lastwarn});
    end
    
    if f2
        warning ('S2_Angles:PortsOrder','Suction close and delivery open are too close. Open system, process unfeasible')
        fOK = 0;
        SX_Logfile ('e',{lastwarn});
    end
    
    
    % check for nozzle position
    for i_nz=1:length(NOZZLES.theta_nz)
                
        if NOZZLES.theta_nz(i_nz) < 0
            warning('S2_Angles:SuctionNzls','Active nozzle (%i): Nozzle angle must be higher or equal to 0', i_nz);
            SX_Logfile ('e',{lastwarn});
            fOK = 0;
        end
        
        if NOZZLES.theta_nz(i_nz) >= th_DisOpen-Gamma/2
            warning('S2_Angles:DeliveryNzls','Active nozzle (%i): SVEC does not support Nozzles during delivery. Nozzle angle must be less than %2.0f°', i_nz, rad2deg(th_DisOpen-Gamma/2));
            SX_Logfile ('e',{lastwarn});
            fOK = 0;
        end       
    end
    
    clear f1 i_nz
end