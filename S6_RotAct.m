function [thv_res,Tx,Ty,Mx,My,Mz,pck1,pck2,pck3,pck4,pck5] = S6_RotAct(xsv,L,d1,bB,R_x,R_xg,R_y,R_yg,qG,C_bronz,C_rot,theta,Npt_cell,shm,zeta,rot_dir,toll_t,n_van,pos_SucOpen,pos_SucClose,pos_DisOpen,pos_DisClose)
% Evaluates Tx, Ty, Mx, My and Mz actions along the shaft. 
% Refer to Lazzari for more information about conventions, referments and
% nomenclature.
%
% INPUT
% xs [m]        : shaft mesh points values [Npt x shm]                                     
% L  [m]        : rotor length                                         
% d1 [m]        : distance from bushing to rotor                     
% bB [m]        : bushing length                                     
% R_x [N]       : overall x-component of shaft reaction force
% R_xg [N]         : x-component of the total resultant on the rotor accounting for the rotor-shaft weight force
% R_y [N]       : overall y-component of shaft reaction force
% R_yg [N]      : y-component of the total resultant on the rotor accounting for the rotor-shaft weight force
% qG  [N/m]     : shaft weigth as distributed load                     
% C_bronz [Nm]  : bushings resistance torque                           
% C_rot   [Nm]  : overall torque acting on shaft due to reaction forces                   
% theta [rad]   : theta angle [-gamma ; 2*pi]                   
% Npt_cell [-]  : number of discretization points between 2 vanes                               
% shm [m]       : shaft length mesh points 
% zeta [°]         : angles of x axis with respect to X
% rot_dir [-]      : direction of rotation     1: clockwise       2: counterclockwise
% toll_t [-]    : theoric tolerance (default: 1e-9)
%
% OUTPUT 
% thv_res [rad] : theta values after resampling   [1 x Npt*]
% Tx [N]        : Tx action evolution in theta    [Npt x shm]
% Ty [N]        : Ty action evolution in theta    [Npt x shm]
% Mx [Nm]       : Mx action evolution in theta    [Npt x shm]
% My [Nm]       : My action evolution in theta    [Npt x shm]
% Mz [Nm]       : Mz action evolution in theta    [Npt x shm]
% pck1-pck5 [-] : position checks of the mesh points along the shaft   [Npt x shm]
% 
% NOTE:
% *: Changes its dimension after resampling. Resampling level is not known
%    a priori.

%%  PRELIMINARY CALCULATION %%
% It will be processed a re-sampling on a lower theta mesh point number, in 
% order to keep code performances high.

thv  = theta(Npt_cell+1:end);                                       % theta angle [0 ; 2*pi]
idx_res = S6_Reschos(length(thv),Npt_cell,n_van,pos_SucOpen,pos_SucClose,pos_DisOpen,pos_DisClose);       % Discretization step along the rotor

thv_res = thv(idx_res);                      % theta values (re-sampling).
R_x_res = R_x(idx_res);                      % pressure and blade forces, x component (re-sampling).
R_y_res = R_y(idx_res);                      % pressure and blade forces, y component (re-sampling).
RB1_x   = -R_xg(idx_res)/2;                  % load on 1st bushing, x component (re-sampling).
RB1_y   = R_yg(idx_res)/2;                   % load on 1st bushing, y component (re-sampling).
RB2_x   = RB1_x;                             % load on 2nd bushing, x component (re-sampling).
RB2_y   = RB1_y;                             % load on 2nd bushing, y component (re-sampling).
CB1     = C_bronz(idx_res)/2;                % 1st bushing friction torque (re-sampling).
CB2     = CB1;                               % 2nd bushing friction torque (re-sampling).
Cm1     = C_rot(idx_res);                    % ideal mechanical torque (re-sampling).

xs      = (repmat(xsv',1,length(thv_res)))';    % shaft mesh grid.
R_xs    = repmat(R_x_res',1,shm);               % pressure and blade forces, x component (filling).
R_ys    = repmat(R_y_res',1,shm);               % pressure and blade forces, y component (filling).
RB1_xs  = repmat(RB1_x',1,shm);                 % load on 1st bushing, x component (filling).
RB1_ys  = repmat(RB1_y',1,shm);                 % load on 1st bushing, y component (filling).
RB2_xs  = repmat(RB2_x',1,shm);                 % load on 2nd bushing, x component (filling).
RB2_ys  = repmat(RB2_y',1,shm);                 % load on 2nd bushing, y component (filling).
CB1s    = repmat(CB1',1,shm);                   % 1st bushing friction torque (filling).
CB2s    = repmat(CB2',1,shm);                   % 2nd bushing friction torque (filling).
Cm1s    = repmat(Cm1',1,shm);                   % ideal mechanical torque (filling).

qGs     = repmat(qG',1,length(Cm1));            % rotor weight filling.
qGs     = qGs';                                 % rotor weight format correction.

clear R_x_res R_y_res RB1_x RB1_y RB2_x RB2_y

% "pck" logical matrices define the position assumed by mesh points along the
% shaft. They return 1 if the considered mash point is in the specificated
% range of positions (0 otherwise).

pck1 = (xs <= bB) + 0;
pck2 = (xs <= (bB + d1)) - pck1;
pck3 = (xs <= (bB + d1 + L)) - pck1 - pck2;
pck4 = (xs <= (bB + d1 + L + d1)) - pck1 - pck2 - pck3;
pck5 = 1 - pck1 - pck2 - pck3 - pck4;

% terms that often appear in formulas.

trm1 = (bB/2)*RB1_xs;
trm2 = (bB/2)*RB1_ys;
trm3 = R_xs;
trm4 = R_ys;
trm5 = trm3*(bB + d1 + (L/2));
trm6 = trm4*(bB + d1 + (L/2));
lgt1 = xs - (bB + d1);
lgt2 = xs - (bB + d1 + L + d1);
lgt3 = xs + (bB + d1);
lgt4 = xs + (bB + d1 + L + d1);

%% INNER ACTIONS EVALUATION %%

Tx1 = ((1/bB)*RB1_xs).*xs;
Tx2 = RB1_xs;
Tx3 = RB1_xs + ((1/L)*trm3).*lgt1;
Tx4 = RB1_xs + trm3;
Tx5 = Tx4 + ((1/bB)*RB2_xs).*lgt2;
Tx  = Tx1.*pck1 + Tx2.*pck2 + Tx3.*pck3 + Tx4.*pck4 + Tx5.*pck5 - qGs*sin(zeta).*xs;

Ty1 = ((1/bB)*RB1_ys).*xs;
Ty2 = RB1_ys;
Ty3 = RB1_ys - ((1/L)*trm4).*lgt1;
Ty4 = RB1_ys - trm4;
Ty5 = Ty4 + ((1/bB)*RB2_ys).*lgt2;
Ty  = Ty1.*pck1 + Ty2.*pck2 + Ty3.*pck3 + Ty4.*pck4 + Ty5.*pck5 + qGs*rot_dir*cos(zeta).*xs;

My1 = (((1/bB)*RB1_xs).*((1/2)*xs) - Tx).*xs;
My2 = trm1 - (Tx.*xs);
My3 = trm1 + (1/L)*trm3.*lgt1.*((1/2)*lgt3) - (Tx.*xs);
My4 = trm1 + trm5 - (Tx.*xs);
My5 = trm1 + trm5 + ((1/bB)*RB2_xs).*lgt2.*((1/2)*lgt4) - (Tx.*xs);
My  = My1.*pck1 + My2.*pck2 + My3.*pck3 + My4.*pck4 + My5.*pck5 - qGs*rot_dir*sin(zeta).*(xs.^2)./2;

Mx1 = (((1/bB)*RB1_ys).*((1/2)*xs) - Ty).*xs;
Mx2 = trm2 - (Ty.*xs);
Mx3 = trm2 - (1/L)*trm4.*lgt1.*((1/2)*lgt3) - (Ty.*xs);
Mx4 = trm2 - trm6 - (Ty.*xs);
Mx5 = trm2 - trm6 + ((1/bB)*RB2_ys).*lgt2.*((1/2)*lgt4) - (Ty.*xs);
Mx  = Mx1.*pck1 + Mx2.*pck2 + Mx3.*pck3 + Mx4.*pck4 + Mx5.*pck5 + qGs*cos(zeta).*(xs.^2)./2;

Mz1 = ((1/bB)*CB1s).*xs;
Mz2 = CB1s;
Mz3 = CB1s + lgt1.*((1/L)*Cm1s);
Mz4 = CB1s + Cm1s;
Mz5 = (CB1s + Cm1s) + ((1/bB)*CB2s).*lgt2;
Mz  = Mz1.*pck1 + Mz2.*pck2 + Mz3.*pck3 + Mz4.*pck4 + Mz5.*pck5;

clear trm1 trm2 trm3 trm4 trm5 trm6 lgt1 lgt2 lgt3 lgt4
%% INTERNAL CHECK %%
% check on the correct value of internal actions at the terminal section of
% the shaft

err_ra1 = Tx(:,shm);
err_ra2 = Ty(:,shm);
err_ra3 = Mx(:,shm);
err_ra4 = My(:,shm);
err_ra5 = Mz(:,shm) - (Cm1 + CB1 + CB2)';

chck = max([max(abs(err_ra1)) max(abs(err_ra2)) max(abs(err_ra3)) max(abs(err_ra4)) max(abs(err_ra5))]);

if chck > toll_t
    warning('S6_RotAct:inner_actions','Rotor inner action congruence is not verified.')
    SX_Logfile ('d',{lastwarn});
end

end

%% SUBFUNCTION S6_Reschos
% Choses a reasonable level of resampling (RESampling CHOSer).

function [idx] = S6_Reschos(Npt,Npt_cell,n_van,pos_SucOpen,pos_SucClose,pos_DisOpen,pos_DisClose)
  % ====================================================================
  % ====================================================================
  % Reschos takes the number of mesh points (Npt) in input and returns a
  % vector of indexes accounting both for the resampling and for the most
  % critical angular positions. The latter represent the positions of the
  % refernce 'following' blade where at the same time another blade is
  % passing over the discharge opening. This is due to the fact that in the
  % last vane the pressure is lower than the one at the discharge opening
  % increasing locally the stress value. To better understand this, please
  % refer to FRANZETTI-PERSICO, analysisng the periodic trend of the
  % equivalent stress.
  %
  % INPUT:
  % Npt          [-]  : number of grid discretization points in [0:2pi]
  % Npt_cell     [-]  : number of discretization points between 2 vanes
  % n_van        [-]  : number of vanes
  % pos_DisOpen  [-]  : discretized delivery open angle position in array theta
  % pos_DisClose [-]  : discretized delivery close angle position in array theta
  % process      [-]  : logic parameter (1 compression / 2 expansion)
  %
  % OUTPUT
  % idx          [-]  : indexes of discretization points
  
     fctz = factor(Npt);
     fctz(fctz >= 150) = [];
         
     val = cumprod(fctz)*max(fctz);
     val = min(val(val >= 100 & val <= 150));
          
     idx_val = 1:val:Npt;
     
     pos_cr1 = (pos_SucOpen-Npt_cell):-Npt_cell:(pos_SucOpen-n_van*Npt_cell);
     pos_cr2 = (pos_SucClose-Npt_cell):-Npt_cell:(pos_SucClose-n_van*Npt_cell);
     pos_cr3 = (pos_DisOpen-n_van*Npt_cell):Npt_cell:(pos_DisOpen-Npt_cell);
     pos_cr4 = (pos_DisClose-n_van*Npt_cell):Npt_cell:(pos_DisClose-Npt_cell);
     idx_cr = [pos_cr1, pos_cr2, pos_cr3, pos_cr4];
     
     idx_cr(idx_cr < 0) = idx_cr(idx_cr < 0) + Npt;    % sign correction accounting for the periodicity of the phenomenon
     idx = unique([idx_val,idx_cr]);                  % sort the indexes vector
          
  % ====================================================================
  % ====================================================================
  
end