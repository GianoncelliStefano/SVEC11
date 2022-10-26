function [sigmaI_max,Fs_traz,Fs_FA] = S6_StateOfStress(MaxMf,N_sMaxMf,T_sMaxMf,s_1D,L,sigmaR_traz,sigmaR_comp,ratio_FA,b2,b3)
% This function computes a stress analysis for the point where the maximum
% bending moment has been found. The analysis is performed by calculatingthe Mohr
% circle and the tensile safety coefficinet for both static and fatigue stress
% INPUT
% MaxMf [Nm]        : absolute maximum of Mf
% N_sMaxMf [N]      : internal action N for the angular position and vane abscissa where maxMf has been found
% T_sMaxMf [N]      : internal action T for the angular position and vane abscissa where maxMf has been found
% s_1D [m]          : vane thickness for 1D geometry
% L [m]             : vane length
% sigmaR_traz [Pa]  : ultimate tensile stress
% sigmaR_comp [pa]  : ultimate compressive stress
% ratio_FA [-]      : fatigue stress ratio
% b2 [-]            : fatigue size factor
% b3 [-]            : fatigue surface finish factor
%
% OUTPUT
% sigmaI_max [pa]   : greatest I main stress obtained via Mohr circle
% Fs_traz [-]       : static tensile safety factor
% Fs_FA [-]         : fatigue safety factor
% 
% NOTE:  Since, cast iron is a brittle material, fatigue test is performed
%         only for tensile state of stress. The test is done for the y
%         coordinate corresponding to the maximum sigmaI (posI).
%         If also a shear stress is present where sigmaI is maximum, a mixed
%         state of state occours. SVEC cannot perform such fatigue test.
    %% STATIC TEST %%
    % definition of y coordinate (coordinate across vane thickness)
    y = linspace(-s_1D/2,s_1D/2,10);

    % calculation of strains as function on y coordinates [pa]
    sigma_fl = MaxMf * y/(L*s_1D^3/12);                  % bending
    sigma_N  = N_sMaxMf/(s_1D*L);                        % compression
    tau_xy   = 6*T_sMaxMf*(s_1D^2/4 - y.^2)/(L*s_1D^3);  % shear on xy plane
    sigma_x  = sigma_fl+sigma_N;                         % compression/bending equivalent stress
    clear y sigma_fl sigma_N

    % Mohr circle calculation [pa]
    sigma_o  = (sigma_x)/2;
    R_Mohr   = sqrt(((sigma_x)/2).^2 + tau_xy.^2);  % Mohr circle radius
    sigmaI   = sigma_o + R_Mohr;                    % I main stress as function on y coordinates (biggest main stress)
    sigmaIII = sigma_o - R_Mohr;                    % III main stress as function on y coordinates (smallest main stress)
    
    [sigmaI_max, posI] = max(sigmaI);             % greatest I main stress obtained via Mohr circle
    sigmaIII_min       = min(sigmaIII);           % smallest III main stress obtained via Mohr circle
    clear sigma_x sigma_o R_Mohr sigmaI sigmaIII

    % Static test: Galileo-Rankine criterion of resistance
    Fs_traz  =  sigmaR_traz/sigmaI_max;    % static tensile safety factor
    Fs_comp  = -sigmaR_comp/sigmaIII_min;  % static compression safety factor

    %% FATIGUE TEST %%
    % fatigue test performed only if no shear stress is present (see NOTE)
    if tau_xy(posI) == 0 
        % fatigue cycle parameters [pa]
        sigma_max  = sigmaI_max;
        sigma_min  = 0;
        sigma_medP = (sigma_max+sigma_min)/2;
        sigma_altP = (sigma_max-sigma_min)/2;
        
        % Haigh diagram for fatigue stress assesment: the coordinates of Q
        % point are found as interseption between limit curve and operating
        % curve. II order equation.
        sigmaFA = ratio_FA *sigmaR_traz*b2*b3;   % ultimate fatigue stress if sigma_medP is null
        % II order equation coefficient
        a = 1;
        b = sigmaR_traz +sigmaFA;
        c = -sigmaFA *sigmaR_traz;
        delta = b^2 -4*a*c;
        sigma_medQ = (-b +sqrt(delta))/(2*a);
        clear a b c delta sigmaFA sigmaFA sigma_max sigma_min 

        % OP and OQ segment of haigh diagram are calculated
        OP = sqrt(sigma_medP^2 + sigma_altP^2);
        OQ = sqrt(sigma_medQ^2 + sigma_medQ^2);
        
       % Fatigue safety coefficnet is calculated
        Fs_FA = OQ/OP;
        clear sigma_medP sigma_altP sigma_medQ OP OQ

    else
        warning('S5_StressTest:fatigue','Shear stress is present. It is necessary a Fatigue test for mixed state of stress')
        SX_Logfile ('w',{lastwarn});
    end

end
