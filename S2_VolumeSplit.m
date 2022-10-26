function [V_suc,V_dis,V_comp] = S2_VolumeSplit (V_cell,pos_SucClose,pos_DisOpen,Npt_cell)
% This function split the volume array in 3  relevant sub array
% INPUT
% V_cell [m3]       : cell volumes array of a complete vane revolution (-gamma:2pi)
% pos_SucClose [-]  : discretized suction close angle position in array theta
% pos_DisOpen [-]   : discretized delivery open angle position in array theta
% Npt_cell [-]      : number of discretization points between 2 vanes
%
% OUTPUT
% V_suc [m3]   : volume array of suction phase
% V_dis [m3]   : volume array of discharge phase
% V_comp [m3]  : volume array of closed chamber phase

    %% COMPUTATIONS %%
    V_suc   = V_cell(1:pos_SucClose-1);
    V_comp  = V_cell(pos_SucClose:pos_DisOpen-Npt_cell);
    V_dis   = V_cell(pos_DisOpen-Npt_cell+1:end-1);
    

end