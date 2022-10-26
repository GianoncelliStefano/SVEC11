function [p_out,p_cell,p_1,p_2,p_geom] = S3_Pressure(theta,theta_vane,p_comp,pos_SucClose,pos_DisOpen,Npt_cell,p_del)
% It determines the trend of pressures in the reference volume for a complete rotation [0:2pi]
% INPUT 
% theta [rad]       : angular position vector of reference vane [-gamma:2pi]
% theta_vane [rad]  : angular position vector [0:2pi]
% p_comp [Pa]       : vector of pressures from suction closure to discharge aperture
% pos_in [-]        : position of the suction closure angle in vector theta                 
% pos_out [-]       : position of the discharge aperture angle in vector theta                 
% Npt_cell [-]      : number of points of discretization of a cell
%
% OUTPUT
% p_van [Pa]  : pressure trend for a complete turning of rotor [0:2pi]         
% p_out [Pa]  : pressure at end of pre-compression 
% p_cell [Pa] : pressure vector of the cell [-gamma:2pi]
% p_2 [pa]    : high pressure acting on vane in [0:2pi]
% p_1 [pa]    : low pressure acting on vane in [0:2pi]
% p_geom [pa] : geometrical pressure (pressure in the cell before delivery opening)
%
% NOTE: 
% 1) the reference van for p_van in the trailing vane. The pressure
%     correspond to the one on pressure side (compression logic)
% 2) The reference cell for p_cell is the one downstream the trailing vane 
% 3) The pressure of precompression p_out usually does not coincide with
%    the delivery grid one p_del

    %% PRELIMINARY CALCULATIONS %%
    p_van  = NaN(size(theta_vane));  % preallocation of p_van vector
    p_cell = NaN(size(theta));       % Preallocation of p_cell vector [pa]
    
    % it is necessary to translate the ports angle by  gamma (= Npt_cell points)
    % because ports positions are computed from theta=-gamma, while p_comp is computed from theta=0
    pos_in_p  = pos_SucClose - Npt_cell;  % position of the suction closure angle in vector theta_vane 
    pos_out_p = pos_DisOpen - Npt_cell;   % position of the discharge aperture angle in vector theta_vane

    %% PRESSURE ON VANE %%
    p_geom = p_comp(end);                          % pressure from geometric process
    
    % pressure vector is created [pa]
    p_van(1:pos_in_p-1) = p_comp(1);               % Suction -  theta = 0 : suction close 
    p_van(pos_in_p:pos_out_p-Npt_cell) = p_comp;   % Compression - theta = suction close : delivery open
    if isnan(p_del)                                % Delivery - theta = delivery close : 2pi
        p_van(pos_out_p-Npt_cell+1:end) = p_geom;  % pressure from geometric process
    else
        p_van(pos_out_p-Npt_cell+1:end) = p_del;   % pressure from user input
    end

    % correction on last element
    p_van(end) = p_van(1);  % close PV graph
    p_out = p_van(end-1);   % end compression pressure (is either p_geom or p_del)
    
    %% CELL PRESSURE %% 
    p_cell(1:Npt_cell)     = p_van(1); % suction side [-gamma:0]
    p_cell(Npt_cell+1:end) = p_van;    % compresion [0:2pi]
    p_cell(end)            = p_van(1); % to close p-V diagram
    
    %% PRESSURE ON PRESSURE AND SUCTION SIDE OF VANE
    p_2             = p_van;                   % pressure side
    p_1             = circshift(p_2,Npt_cell); % suction side (needs correction)
    p_1(1:Npt_cell) = p_2(1:Npt_cell);         % resetting first values
    p_2(end)        = p_1(end);                % to close PV diagram

end