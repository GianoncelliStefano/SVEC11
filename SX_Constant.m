function [varargout] = SX_Constant(K)
% Engineering constants in S.I units
%       CONSTANT      |    UNITS
% ----------------------------------
%      UniGasConstant |  [kJ/kmol/K]
%      Avogadro       |  [mol^-1] 
%      Gravity        |  [m/s^2]
%      Boltzman       |  [W/m^2/K^4]
%      LightSpeed     |  [m/s]
% INPUT
%   K: p-cell array with required constants i.e {'const1','const2',...}
% OUTPUT
%   (1-by-p) array with retrieved values    i.e [value1, value2, value... ]
% EXAMPLE
%   [R,k]=Constant({'UniGasConstant','Avogadro'})
% HISTORY
%   September 2017 - Lopez Juan P.

    varargout = cell(size(K));
    for j=1:length(K)
        switch K{j}
            case 'UniGasConstant'
                varargout{j} = 8.314472 ;     % [J/mol/K]
            case 'Avogadro'
                varargout{j} = 6.022141e23;   % [mol^-1]
            case 'Gravity'
                varargout{j} = 9.81;          % [m/s^2]
            case 'Boltzman'
                varargout{j} = 5.670367e-8;   % [W/m^2.K^4]
            case 'LightSpeed'
                varargout{j} = 2.99792458e8;  % [m/s]
        end
    end
end