function [sVB,sBI,pos_contct] = S6_SaxisVane(BI,VB,stp,Npt)
% This function performs a discretization of the vane lenght for each angular position
% INPUT
% BI [m]          : excursion of 1d vane (corresponds to BI segment)
% VB [m]          : segment from bottom vane to vane-rotor interseption
% stp [m]         : vane discretization step for stress analysis
% Npt [-]         : number of grid discretization points in [0:2pi]
%
% OUTPUT
% sVB [-]         : cell-array of DB discretized lenght 
% sBI [-]         : cell-array of BC discretized lenght
% pos_contct [-]  : position of the point preceeding rotor - vane contact point for each discretized angular position
%
% NOTE:  - For each discretized angular position, vane abscissa s is computed.
%        - Each element of the cell array has different dimension.
%        - The dimension of each element is given by the vane discretization step,
%           and is rounded up to the next integer

    %% CALCULATIONS %%
    % vane dimension for discretization are computed
    BI    = max(BI,0);  % negative values are set to zero

    % Preallocation of cell arrays
    sVB        = cell(1,Npt);
    sBI        = cell(1,Npt);
    pos_contct = NaN(1,Npt);

    % construction of vane abscissa
    for i=1:Npt
        sVB{1,i} = linspace(0,VB(i),ceil(VB(i)./stp));
        sBI{1,i} = linspace(0,BI(i),ceil(BI(i)./stp));
        
        % Because of discretization, the rotor-vane contact point (point B) is splitted in two
        pos_contct(1,i) = length(sVB{1,i});   % point preceeding rotor-vane contact position
        
        % empty values are set to zero
        if isempty(sVB{1,i})
            sVB{1,i} = 0;
        end
        if isempty(sBI{1,i})
            sBI{1,i} = 0;
        end
    end

end