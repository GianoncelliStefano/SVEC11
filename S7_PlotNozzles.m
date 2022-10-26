function S7_PlotNozzles(T_l,T_lmean,V_inj_cls,Dg,th_in,th_out,Gamma,NOZZLES)
% This function creates the plots of nozzles oil injected and oil
% temperature during the closed cell phase
% INPUT
% T_l [K]        : liquid temperature during closed chamber phase
% T_lmean [K]    : oil adiabatic mixing temperature for each discretization step
% V_inj_cls [m3] : volume of oil injected by each nozzle class
% Dg [m]         : droplet diameters for each class
% th_in [rad]    : suction close angle (takes into account if user wants to maximise suction volume & expansion ports flipping)
% th_out [rad]   : delivery open angle (takes into account expansion ports flipping)
% Gamma [rad]    : angular extention of a cell (cavity between 2 vanes)
% NOZZLES        : structure with nozzles parameters

    %% DEFINITIONS %%
    fOIL  = 0;           % default: no oil is injected 
    
    if ~isempty(T_l)     % oil present
        fOIL      = 1;
        NptPlot   = size(T_l,1);                                      % number of elements on X axis
        x_rad     = linspace(th_in,th_out-Gamma,NptPlot);             % X axis in radiants
%         colmatrix = SX_DistinguishableColors(size(T_l,2),{'w','k'});  % generates a set of complete different colors
% not suitable command with Matlab 2017
        colmatrix = colormap;

    end

    if fOIL        
        %% CONVERSION %%
        x_deg     = rad2deg(x_rad);             % [rad -> deg]
        theta_nz  = rad2deg(NOZZLES.theta_nz);  % [rad -> deg]
        Gamma     = rad2deg(Gamma);             % [rad -> deg]
        V_inj_cls   = V_inj_cls*1e9;                % [m3 -> mm3]
        Dg        = Dg*1e6;                     % [m -> microm]
        T_l       = T_l-273.15;                 % [K -> °C]
        T_lmean   = T_lmean-273.15;             % [K -> °C]
        Tl_in     = NOZZLES.Tl_in-273.15;       % [K -> °C]

        %% PLOT OIL TEMPERATURE %%
        figure('Name','Fig_4 Nozzles');
        subplot(3,1,1);
        set(gca, 'ColorOrder',colmatrix, 'NextPlot', 'replacechildren');   % assign color list
        plot(x_deg,T_lmean,'k--',x_deg,T_l'), hold on;                     % plot oil temperature
        plot((theta_nz-Gamma/2),Tl_in,'dk');                               % plot simulation oil injection angle
        %options
        grid on, box on;
        xlim('auto');
        title('Oil temperature during closed cell phase, for each droplet class');
        xlabel('Trailing vane angle [deg]');
        ylabel('Temperature [°C]');
        legend('T_{MIX adb}','Location','west');

        %% PLOT VOLUME OF OIL INJECTED %%
        subplot(3,1,2);
        [~,Vcls,Xtick,dim] = NozzlesPro (NOZZLES,V_inj_cls);   % create adjusted vector of volumes 
        bVol = bar(Vcls,0.3);                                % plot volume of oil injected
        %options
        AssignColor(NOZZLES,colmatrix,bVol,dim);             % assign different colors to each bar element
        set(gca,'xticklabel',NOZZLES.type_nz);               % set X axis names
        set(gca,'XTick',Xtick);                              % set X axis names positions
        grid on, box on;
        xlim('auto');
        title('Volume of injected oil for each droplet class');
        xlabel('Nozzle classes');
        ylabel('Volume [mm3]');

        %% PLOT DROPLETS DIAMETER %%
        subplot(3,1,3);
        [~,Dcls,Xtick,dim] = NozzlesPro (NOZZLES,Dg);  % create adjusted vector of diameters
        bDia = bar(Dcls,0.3);                          % plot diameter of droplets of oil injected
        %options
        AssignColor(NOZZLES,colmatrix,bDia,dim)        % assign different colors to each bar element
        set(gca,'xticklabel',NOZZLES.type_nz);         % set X axis names
        set(gca,'XTick',Xtick);                        % set X axis names positions
        grid on, box on;
        xlim('auto');
        title('Diameter of  each droplets class');
        xlabel('Nozzle classes');
        ylabel('Diameter [\mum]');

    end

end


function [NameOUT,VarOUT,Xtick,dim] = NozzlesPro (NOZZLES,VarIN)
% This function creates a vector containing the VarIN data separated by a
% NaN, so that each group of data correspond to a single nozzle
% in this way is possible to plot the bar chart, 'cheating' Matlab,
% and obtain a kind of group separation for each class. The issue is that
% the x-axis will not be grouped
% INPUT
% NOZZLES   : structure with nozzles parameters
% VarIN     : row array with a property of the nozzles
%
% OUTPUT
% NameOUT  : cell array needed for the creation of the bar chart
% VarOUT   : nozzles data separated in nozzles groups by NaN
% Xtick    : Xtick label index, so that it is centered with the bar chart groups
% dim      : length of VarOUT
% NOTE: Even though this function is a mess, it is possible to set the color of each
%       individual bar later with the function AssignColor.
%       To create the 2 bar plots with this function (example for volume bar chart)
%         [Typecls,Vcls,dim] = NozzlesPro (NOZZLES,V_inj_cls);
%         bVol = bar(Vcls,0.5);
%         set(gca,'xticklabel',Typecls);
%         AssignColor(NOZZLES,refplot,bDia,dim)

dim     = length(VarIN)+length(NOZZLES.f_nz)-1;
VarOUT  = NaN(1,dim);
NameOUT = cell(1,dim);
Xtick   = NaN(1,length(NOZZLES.f_nz));
n_nzl   = 1;              % current nozzle
j       = 1;              % index of VarOUT
skip    = 0;              % number of NaN inserted
k       = 1;              % index of class of current nozzle
while j <= dim
    VarOUT(j)  = VarIN(j-skip);
    NameOUT{j} = NOZZLES.type_nz{n_nzl};
    if k == NOZZLES.num_cls(n_nzl) && j ~= dim
        Xtick(n_nzl) = j-k/2+0.5;
        VarOUT(j+1)  = NaN;
        NameOUT{j+1} = '*';
        n_nzl        = n_nzl + 1;
        j            = j+2;
        skip         = skip + 1;
        k = 1;
    else
        j = j+1;
        k = k+1;
    end
end
Xtick(end) = (j-1)-(k-1)/2+0.5;
end


function AssignColor (NOZZLES,colmatrix,plotIN,dim)
% This function can copy the color in the RGB triplet matrix colmatrix and assgn
% them to the bar chart PLOTIN. This function must be used with NozzlePro
% INPUT
% NOZZLES   : structure with nozzles parameters
% colmatrix : matrix of RGB triplets color (e.g. [0.2, 0.65, 0.8; 1, 0.4, 0.4])
% plotIN    : bar chart to assign the color from refplot
% dim       : numer of element on Plot IN (some of them are NaN)
%
% NOTE: Even though this function is a mess, it is the only way to set the color of each
%       individual bar created with NozzlesPro, so that the color matches
%       with the one of the temperature plot.
%       To create the 2 bar plots with this function (example for volume bar chart)
%         [Typecls,Vcls,dim] = NozzlesPro (NOZZLES,V_inj_cls);
%         bVol = bar(Vcls,0.5);
%         AssignColor(NOZZLES,colmatrix,bDia,dim)
%         set(gca,'xticklabel',Typecls);
%         set(gca,'XTick',Xtick);  {alternatevely is possible to use: set(gca,'XTick',NOZZLES.type_nz)}

plotIN.FaceColor = 'flat';
n_nzl   = 1;             % current nozzle
j       = 1;             % index of line elements in refplot
skip    = 0;             % number of NaN skipped
k       = 1;             % index of class of current nozzle
while j <= dim
    plotIN.CData(j,:) = colmatrix((j-skip),:);
    if k == NOZZLES.num_cls(n_nzl) && j ~= dim
        n_nzl       = n_nzl + 1;
        j           = j+2;
        skip        = skip + 1;
        k           = 1;
    else
        j = j+1;
        k = k+1;
    end
end
end


