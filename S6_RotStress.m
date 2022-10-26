function [sigmaeqD,rSFsD,rSFfD,xscrD,thcrD,sigmaeqDsh,rSFsDsh,rSFfDsh,xscrDsh,thcrDsh] = S6_RotStress(thv_res,xsv,pck1,pck2,pck3,pck4,pck5,d_hub,d,l,Tx,Ty,Mx,My,Mz,Rmt_s,Rp02_s,sigma_fas,b2_s,b3_s)
% Rotor static and fatigue stress analysis. See Lazzari for more
% informations about conventions, referments and nomenclature.
%
% INPUT 
% thv_res    [rad] : theta values after resampling                       
% xsv        [m]   : rotor+shaft axial length discretization                            
% pck1-pck5  [-]   : position checks of the mesh points along the shaft         
% d_hub      [m]   : shaft diameter                                      
% d          [m]   : rotor diameter                                      
% l          [m]   : blade height                                         
% Tx         [N]   : Tx action evolution in theta                       
% Ty         [N]   : Ty action evolution in theta                        
% Mx         [Nm]  : Mx action evolution in theta                       
% My         [Nm]  : My action evolution in theta                       
% Mz         [Nm]  : Mz action evolution in theta                       
% Rmt_s      [MPa] : Ultimate tensile strength                    
% Rp02_s     [MPa] : Yield strength                                
% sigma_fas  [MPa] : Fatigue strength 
% b2_s       [-]   : Dimension factor 
% b3_s       [-]   : Surface roughness factor 
%
% OUTPUT 
% sigmaeqD   [MPa] : absolute maximum  shaft+rotor stress value (ductile analysis)    
% rSFsD      [-]   : shaft+rotor ductile static safety factor                 
% rSFfD      [-]   : shaft+rotor shaft ductile fatigue safety factor              
% xscrD      [m]   : critical position on the shaft+rotor (ductile analysis)   
% thcrD      [rad] : critical shaft+rotor position in theta (ductile analysis) 
% sigmaeqDsh [MPa] : absolute maximum  shaft stress value (ductile analysis)  
% rSFsDsh    [-]   : shaft ductile static safety factor  
% rSFfDsh    [-]   : shaft shaft ductile fatigue safety factor 
% xscrDsh    [m]   : critical position on the shaft (ductile analysis)   
% thcrDsh    [rad] : critical shaft position in theta (ductile analysis) 
%
% NOTE: Stress matrix was firstly designed by D. Coletta (see Shaft Stress code)
%       according to Lazzari thesis, this version can be found in the last lines. 
%       The new stress matrix has been defined by Franzetti-Perisco (2019).

%% USEFUL DEFINITIONS %%
% Definition of several units which often appear in formulas.
d_trend = d_hub*pck1 + d_hub*pck2 + (d - 2*l)*pck3 + d_hub*pck4 + d_hub*pck5; % Section diameter trend.
r = d_trend./2;     % Section radius trend.
J = (pi/4)*(r.^4);  % Inertia moment.
S = pi*r.^2;        % Section trend.

% Mesh definition.
[thm,shm]   = size(Tx);

%% STATIC ANALYSIS %%
% Pre-allocation.
sigmamaxD   = zeros(thm,shm); % Ductile equivalent stresses.

% In every point theta and in every point along the shaft it must be
% determined the stress matrix "sigma". Stress components in the principal
% reference can be evaluated from the eigenvalues of "sigma". These are
% collected for every point in theta and for every point along the shaft
% in the apposites matrices.

% Ductile analysis for each point along the shaft...    
for i = 1:shm
    % ...and in each point in theta.
    for j = 1:thm
        % Definition of the stress matrix.
        [sigma] = [0                                                   0                                                      4/3*Tx(j,i)/S(j,i)+0.5*Mz(j,i)*r(j,i)/J(j,i);
                   0                                                   0                                                      4/3*Ty(j,i)/S(j,i)+0.5*Mz(j,i)*r(j,i)/J(j,i);
                   4/3*Tx(j,i)/S(j,i)+0.5*Mz(j,i)*r(j,i)/J(j,i)        4/3*Ty(j,i)/S(j,i)+0.5*Mz(j,i)*r(j,i)/J(j,i)           sqrt((Mx(j,i))^2+(My(j,i))^2)*r(j,i)/J(j,i)];
        
        % Conversion of components in the main referment.
        sigma   = sigma./1000000;   % Conversion from Pa to MPa.
        eigv    = eig(sigma);       % Eigenvalues extraction.
        sigma1  = max(eigv);        % First component in main referment.
        sigma2  = eigv(2);			% Second component in main referment.
        sigma3  = min(eigv);        % Last component in main referment.
        % Equivalent stresses.
        sigmamaxD(j,i) = (1/sqrt(2))*sqrt((sigma2 - sigma1).^2 + (sigma3 - sigma1).^2 + (sigma3 - sigma2).^2);   % Von Mises equivalent stress.
    end
    
end

% Critical values extraction (for the overall rotor and shaft assembly).
% Ductile case.
[val,pos]          = max(sigmamaxD);  % "val" contains the maximum values of each column, "pos" contains their positions.
[sigmaeqD,indxD]   = max(val);        % Extracting the absolute maximum stress value and its position along the shaft.
indthD             = pos(indxD);      % Extracting the position in theta.
xscrD = xsv(indxD);                   % Conversion from index to length.
thcrD = thv_res(indthD);              % Conversion from index to angle.
rSFsD = Rp02_s/sigmaeqD/10^6;         % Ductile static safety factor.

% Critical values extraction (for shaft section only).
% Ductile case.
sigmamaxDsh          = sigmamaxD.*(pck1+pck2+pck4+pck5);
[val,pos]            = max(sigmamaxDsh);  % "val" contains the maximum values of each column, "pos" contains their positions.
[sigmaeqDsh,indxDsh] = max(val);          % Extracting the absolute maximum stress value and its position along the shaft.
indthDsh             = pos(indxDsh);      % Extracting the position in theta.
xscrDsh = xsv(indxDsh);                   % Conversion from index to length.
thcrDsh = thv_res(indthDsh);              % Conversion from index to angle.
rSFsDsh = Rp02_s/sigmaeqDsh/10^6;         % Ductile static safety factor.


%% FATIGUE ANALYSIS
% Start processing (for rotor and shaft assembly).

lhD   = NaN*ones(thm,3);                     % Ductile load history matrix pre-allocation.
%thlhD    = wrapTo2Pi(thv_res - thv_res(indthD));
thlhD = mod(thv_res - thv_res(indthD),2*pi); % Ductile load history angle. Its value is zero whereas theta is equal to "thcrD".

% Ductile analysis for each point in theta
    for j = 1:thm
        % Definition of the rotative bending stress matrix in the most stressed section.
        [sigma] = [0                                                                                0                                                                               cos(thlhD(j))*4/3*Tx(j,indxD)/S(j,indxD)+0.5*Mz(j,indxD)*r(j,indxD)/J(j,indxD) ;
                   0                                                                                0                                                                               cos(thlhD(j))*4/3*Ty(j,indxD)/S(j,indxD)+0.5*Mz(j,indxD)*r(j,indxD)/J(j,indxD) ;
                   cos(thlhD(j))*4/3*Tx(j,indxD)/S(j,indxD)+0.5*Mz(j,indxD)*r(j,indxD)/J(j,indxD)   cos(thlhD(j))*4/3*Ty(j,indxD)/S(j,indxD)+0.5*Mz(j,indxD)*r(j,indxD)/J(j,indxD)  sqrt((cos(thlhD(j))*Mx(j,indxD))^2+(cos(thlhD(j))*My(j,indxD))^2)*r(j,indxD)/J(j,indxD)];          
        sigma   = sigma./1000000;              % Conversion from Pa to MPa.
        lhD(j,:) = sort(eig(sigma),'descend'); % Load history matrix.
    end

% Critical values extraction.
% Both cases.
sigmaeqaD  = (1/sqrt(2))*(sqrt((lhD(:,1) - lhD(:,2)).^2 + (lhD(:,1) - lhD(:,3)).^2 + (lhD(:,2) - lhD(:,3)).^2));  % Von Mises equivalent stress.
sigmaeqamD = mean(sigmaeqaD);                            % Mean stress. 
sigmaeqaaD = max(sigmaeqaD) - sigmaeqamD;                % Amplitude of the stress cycle.
sigma_fa1D = sigmaeqaaD/(1-(sigmaeqamD/(Rmt_s/10^6))^2); % Gerber parabola equation 
rSFfD      = sigma_fas/sigma_fa1D/10^6;                  % Ductile fatigue safety factor. 

% Start processing (for shaft only).

lhDsh   = NaN*ones(thm,3);                       % Ductile load history matrix pre-allocation.
%thlhD    = wrapTo2Pi(thv_res - thv_res(indthD));
thlhDsh = mod(thv_res - thv_res(indthDsh),2*pi); % Ductile load history angle. Its value is zero whereas theta is equal to "thcrD".

% Here is evaluated the ductile case.
    for j = 1:thm
        % Definition of the rotative bending stress matrix in the most stressed section.
        [sigma] = [0                                                                                           0                                                                                           cos(thlhDsh(j))*4/3*Tx(j,indxDsh)/S(j,indxDsh)+0.5*Mz(j,indxDsh)*r(j,indxDsh)/J(j,indxDsh) ;
                   0                                                                                           0                                                                                           cos(thlhDsh(j))*4/3*Ty(j,indxDsh)/S(j,indxDsh)+0.5*Mz(j,indxDsh)*r(j,indxDsh)/J(j,indxDsh) ;
                   cos(thlhDsh(j))*4/3*Tx(j,indxDsh)/S(j,indxDsh)+0.5*Mz(j,indxDsh)*r(j,indxDsh)/J(j,indxDsh)  cos(thlhDsh(j))*4/3*Ty(j,indxDsh)/S(j,indxDsh)+0.5*Mz(j,indxDsh)*r(j,indxDsh)/J(j,indxDsh)  sqrt((cos(thlhDsh(j))*Mx(j,indxDsh))^2+(cos(thlhDsh(j))*My(j,indxDsh))^2)*r(j,indxDsh)/J(j,indxDsh)];                  
        sigma   = sigma./1000000;              % Conversion from Pa to MPa.
        lhDsh(j,:) = sort(eig(sigma),'descend'); % Load history matrix.
    end

% Critical values extraction.
% Both cases.
sigmaeqaDsh  = (1/sqrt(2))*(sqrt((lhDsh(:,1) - lhDsh(:,2)).^2 + (lhDsh(:,1) - lhDsh(:,3)).^2 + (lhDsh(:,2) - lhDsh(:,3)).^2));  % Von Mises equivalent stress.
sigmaeqamDsh = mean(sigmaeqaDsh);                              % Mean stress.
sigmaeqaaDsh = max(sigmaeqaDsh) - sigmaeqamDsh;                % Amplitude of the stress cycle.
sigma_fa1Dsh = sigmaeqaaDsh/(1-(sigmaeqamDsh/(Rmt_s/10^6))^2); % Gerber parabola
rSFfDsh      = sigma_fas/sigma_fa1Dsh/10^6;                    % Ductile fatigue safety factor.

%% FATIGUE SAFETY FACTOR CORRECTION %%
% Safety factor.
rSFfD     = rSFfD*b2_s*b3_s;
rSFfDsh   = rSFfDsh*b2_s*b3_s;
end

% OLD MATRIX DEFINITIONS:
% The new version assumes according to De St. Venant that sigma_x = sigma_y
% = tau_xy = 0 and that the contribution of torque M_z is added to the
% stress component given by shear forces T_x and T_y;
%
% - Static analysis:
%
% ============================================================================================================
% [sigma] = [0                        2*Mz(j,i)*r(j,i)/J(j,i)     3/2*Tx(j,i)/S(j,i);
%            2*Mz(j,i)*r(j,i)/J(j,i)  0                           3/2*Ty(j,i)/S(j,i) ;
%            3/2*Tx(j,i)/S(j,i)       3/2*Ty(j,i)/S(j,i)          sqrt((Mx(j,i))^2+(My(j,i))^2)*r(j,i)/J(j,i)];
% ============================================================================================================
%
% - Fatigue analysis (for the overall rotor and shaft assembly):
%
% ============================================================================================================
% sigma] = [0                                         2*Mz(j,indxD)*r(j,indxD)/J(j,indxD)       cos(thlhD(j))*3/2*Tx(j,indxD)/S(j,indxD) ;
%           2*Mz(j,indxD)*r(j,indxD)/J(j,indxD)       0                                         cos(thlhD(j))*3/2*Ty(j,indxD)/S(j,indxD) ;
%           cos(thlhD(j))*3/2*Tx(j,indxD)/S(j,indxD)  cos(thlhD(j))*3/2*Ty(j,indxD)/S(j,indxD)  sqrt((cos(thlhD(j))*Mx(j,indxD))^2+(cos(thlhD(j))*My(j,indxD))^2)*r(j,indxD)/J(j,indxD)];
% ============================================================================================================ 
%