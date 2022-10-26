function [Fg,qG] = S4_RotWeight(L,d,d_hub,l,s,d1,bB,n_van,rho_s,xsv)
% Rotor+shaft weight force calculation
% INPUT
% L     [m]      : Rotor length
% d     [m]      : Rotor diameter
% d_hub [m]      : Hub diameter
% l     [m]      : Vane height
% s     [m]      : Vane thickness
% d1    [m]      : Rotor end-bushing distance
% bB    [m]      : Bushing length
% n_van [-]      : Number of vanes
% rho_s [kg/m^3] : Shaft density
% xsv   [m]      : Rotor+shaft axial length discretization
%
% OUTPUT
% Fg    [N]     : Global weight force
% qG    [N/m]   : Global weight force as distributed load

   %% PRELIMINARY CALCULATIONS %%
   r_hub	   = d_hub/2;               % [m] Shaft radius 
   r          = d/2;                    % [m] Rotor radius
   g_const = SX_Constant({'Gravity'});
   
   % ====================================================================
   % ====================================================================
   % Method for evaluating qG
   method = 'old';    % 'old': COLLETTA method    'new': Franzetti-Persico method
   % ====================================================================
   % ====================================================================
   
   %% CALCULATIONS %%
   % Rotor+shaft sections determination
   spal  = n_van*l*s;                                                  % [m^2] Slots area 

   pos1 = (xsv <= bB) + 0;                                             % Points along the 1st bushing (flag)
   pos2 = (xsv <= (bB + d1)) - pos1;                                   % Points between the 1st bushing and the rotor (flag)
   pos3 = (xsv <= (bB + d1 + L)) - pos1 - pos2;                        % Points along the rotor (flag)
   pos4 = (xsv <= (bB + d1 + L + d1)) - pos1 - pos2 - pos3;            % Points between the rotor and the 2nd bushing (flag)
   pos5 = 1 - pos1 - pos2 - pos3 - pos4;                               % Points along the 2nd bushing (flag)
   A    = (pi*r_hub^2)*(pos1+pos2+pos4+pos5) + (pi*r^2 - spal)*pos3;   % Rotor+shaft sections trend
   % ====================================================================
   % ==================================================================== 
   clear r r_hub spal pos1 pos2 pos3 pos4 pos5
   % ====================================================================
   % ====================================================================

   %% WEIGHT FORCE CALCULATION ON THE ROTOR+SHAFT %%
   fG     = g_const*rho_s*xsv(2)*A;       % [N] Infinitesimal weight force
   Fg     = sum(fG);                      % [N] Global weight force

   % Global weight force pre-allocation as distributed load
   qG     = NaN(1,length(xsv));  
   
   switch method
       case 'old'
         for i = 1:length(xsv)
           qG(i) = sum(fG(1,1:i));
         end
 % ====================================================================
 % ====================================================================  
        case 'new'
         qG(1,1:length(xsv)) = cumsum(fG);
   end
 % ====================================================================
 % ==================================================================== 
   
   qG = qG./xsv;      % [N/m] Global weight force as distributed load
   
end