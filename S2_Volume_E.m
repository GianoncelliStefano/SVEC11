function V_cell = S2_Volume_E(d,e,L,theta,Gamma,Npt_cell,pos_tang)
% Calculates the volume between two vanes during a complete rotor revolution
% INPUT
% d [m]         : rotor diameter
% e [-]         : ellipse diameter
% L [m]         : rotor length
% theta [rad]   : angular discretized position [-gamma:2pi]
% Gamma [rad]   : angular distance between two vanes
% Npt_cell [-]  : number of discretization points of a cell (space between 2 vanes)
%
% OUTPUT
% V_cell[m^3]   : volumes vector during a complete rotor revolution
%
% NOTE: - theta describes the trailing vane position (pala a seguire)
%       - the reference vane is the trailing vane
% Hypothesis: - radial vane
%             - zero thickness vane

    %% CELL VOLUME CALCULATIONS %%
    % NB: theta is between (-gamma:pi)
    V_cell = NaN(size(theta)); % preallocation of Vvano
    r     = 0.5*d;             % rotor radius [m]
    a     = r/sqrt(1-e^2);     % semi-major axis [m]

    % antiderivative function of cell volume (funzione integranda del volume)
    fvol = @(x)(a*sqrt(1-e^2)./sqrt(1-(e^2*(sin(x).^2)))).^2;  

    % I Integration: theta is between [-gamma:0]
    % Volume between tangency pt and trailing vane [m^3]
    theta_inf = 1e-16;
    for i = 1:Npt_cell
        theta_sup = theta(i)+Gamma;                                  % upper limit of integral
        Arot      = r^2*(theta_sup-theta_inf);                       % area of rotor circular sector
        V_cell(i)  = 0.5*L*(integral(fvol,theta_inf,theta_sup)-Arot);
    end

    % II Integration: theta is between [0:(pi-gamma)]
    % Volume between two vanes [m^3]
    vano_out = length(theta)-Npt_cell;
    for i = pos_tang:vano_out                                        % note: pos_tang = N_vano+1
        theta_inf = theta(i);                                        % lower limit of integral
        theta_sup = theta(i)+Gamma;                                  % upper limit of integral
        Arot      = r^2*(theta_sup-theta_inf);                       % area of rotor circular sector
        V_cell(i)  = 0.5*L*(integral(fvol,theta_inf,theta_sup)-Arot); 
    end

    % III Integration: theta is between [(pi-gamma):pi]
    % Volume between leading vane and tangency pt  [m^3]
    theta_sup = pi+1e-15; % estremo superiore dell'integrale
    for i = vano_out:length(theta)
        theta_inf = theta(i);                                        % lower limit of integral
        Arot      = r^2*(theta_sup-theta_inf);                       % area of rotor circular sector
        V_cell(i)  = 0.5*L*(integral(fvol,theta_inf,theta_sup)-Arot); 
    end

    clear fvol theta_inf theta_sup Arot vano_out
    
    %% CHECKS %%
    % cell volume must be positive
    check1 = V_cell<0;
    if sum(check1) >= 1
        warning('S2_Volume_C:NegativeVolume','Negative cell volume for some angular position!');
        SX_Logfile ('w',{lastwarn});
    end
end
