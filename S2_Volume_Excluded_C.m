function [ Vol_Excluded ] = S2_Volume_Excluded_C( BI, BN, BM, ME, NL, EI, IL, theta_vane, L, s, r_tip, d, AI_1D, AI, AP_i, AP_j, D,toll_t,fDBG)
% This function calculates the three volumes to be sottracted to the limit
% volumes. Three volumes are stored in the variable Vol_excluded (3 x length(theta))
%
% INPUT
% BI [m]                   : real vane excursion measured on vane axis
% BN, BM [m]               : vane-rotor chord
% EI, IL [m]               : vane-tip chord
% ME, NL, [m]              : gas-exposed vane surfaces (lateral surfaces)
% theta_vane [rad]         : discretized angular positions between (0:2pi)
% L [m]                    : compressor axial length
% s [m]                    : vane thickness
% r_tip [m]                : tip radius
% d [m]                    : rotor diameter 
% D [m]                    : stator diameter
% AI_1D, AI, AP_i, AP_j    : geometry elemets (see S2_geometry_C)
%   
% OUTPUT
% Vol_Excluded [m^3]       : variable that stores the three basic volumes
%                             to be sottracted from the limit volume
% Vol_Excluded(1,:)        : vane volume that stands in the FOLLOWING cell - applies to position theta
% Vol_Excluded(2,:)        : vane volume that stands in the PRECEDING cell - applies to position theta + N_pt
% Vol_Excluded(3,:)        : cell volume that stands between vane tip and stator surface
%
%NOTES
% see D Zanchi thesis for references and pictures (Capitolo 2)

%% DECLARATIONS AND CALCULATIONS %%
    % declarations
    r  = d/2; % [m] stator radius
    R  = D/2; % [m] rotor radius 
     
    % preallocation of geometric quantities of quadrilatelars used to  compute
    % the volume of the vane and of the excluded zone aroud the tip volume
    % row 1 and 2 refers to the 2 halves of the vane row 3 refers to zone around the tip
    a   = NaN (3, length(theta_vane)); % [m] first side of quadrilatelar
    b   = NaN (3, length(theta_vane)); % [m] second side of quadrilateral
    c_1 = NaN (3, length(theta_vane)); % [m] first chord of quadrilateral
    c_2 = NaN (3, length(theta_vane)); % [m] second chord of quadrilateral
    h   = NaN (3, length(theta_vane)); % [m] height of quadrilateral
    r_1 = NaN (3, length(theta_vane)); % [m] radius of first chord c_1
    r_2 = NaN (3, length(theta_vane)); % [m] radius of second chord c_2
    
    % construction of geometrical quantities that define the quadrilateral shape
    a(1,:) = ME;
    a(2,:) = BI;
    a_3    = AI_1D - AI; 
    a_3(BI<0) = 0;
    a(3,:) = a_3;
    clear a_3

    b(1,:) = BI;
    b(2,:) = NL;
    b(3,:) = 0;
    c_1(1,:) = BM;
    c_1(2,:) = BN;
    c_1(3,:) = sqrt((AI-AP_i).^2 + AP_j.^2);
    
    c_2(1,:) = EI;
    c_2(2,:) = IL;
    c_2(3,:) = sqrt((AI_1D-AP_i).^2 + AP_j.^2);
    
    h(1,:) = 0.5*s;
    h(2,:) = 0.5*s;
    h(3,:) = AP_j;
    
    r_1(1,:) = r;
    r_1(2,:) = r;
    r_1(3,:) = r_tip;
    
    r_2(1,:) = r_tip;
    r_2(2,:) = r_tip;
    r_2(3,:) = R;
    
    % geometrical quantities correction for theta positions close to tangency
    % the correction is performed using a proportional factor I1 and I2
    Case_1 = (a<0 & b>0);
    Case_2 = (a>0 & b<0);
    Case_3 = (a<=0 & b<=0);
    I1 = b(Case_1)./(b(Case_1) - a(Case_1));
    I2 = a(Case_2)./(-b(Case_2) + a(Case_2));
    
    h(Case_1) = I1.*0.5*s;
    h(Case_2) = I2.*0.5*s;
    h(Case_3) = 0;
    
    if sum(Case_1(1,:))>0
        c_2(Case_1) = I1.*EI;
        c_2(Case_2) = I2.*IL;
        c_1(Case_1) = I1.*BM;
        c_1(Case_2) = I2.*BN;
    else
        c_2(Case_1) = I1.*IL;
        c_2(Case_2) = I1.*EI;
        c_1(Case_1) = I2.*BN;
        c_1(Case_2) = I2.*BM;
    end
    
    c_2(Case_3) = 0;
    c_1(Case_3) = 0;
    a(a<0)      = 0;
    b(b<0)      = 0;
    
    %excluded volumes calculation
    [ Vol_Excluded ] = Volume_Quadrilateral( a ,b ,c_1 ,c_2 ,r_1 ,r_2 ,h ,L );

    if s < toll_t
        Vol_Excluded(:,:) = 0;              %excluded volume is 0 for 1D geometry
    end

    clear a b c_1 c_2 r_1 r_2 h Case_1 Case_2 Case_3
    
%% CHECKS %%
    % the first two rows of Vol_Excluded must be positive
    chck1 = Vol_Excluded(1:2,:)<0;
    if sum(chck1)>=1
        warning('S2_Volume_Excluded_C:VaneVolume','Vane volumes are negative in some angular position');
        SX_Logfile ('w',{lastwarn});
    end, clear chck1
    
    if fDBG ==1
        if r_tip < toll_t && s < toll_t    % then volume_excluded must go to zero
            chck2 = Vol_Excluded > toll_t;
            if sum(chck2)>1
                warning('S2_Volume_Excluded_C:VaneVolume','Vol_Excluded is greater than zero with zero thickness vane');
                SX_Logfile ('d',{lastwarn});
            end, clear chck2
        end
    end
    
end




function [ Vol ] = Volume_Quadrilateral( a ,b ,c_1 ,c_2 ,r_1 ,r_2 ,h ,L )
%function that calculates the volume of the base quadrilateral element,
%given the geometrical parameters previously calculated
%INPUT 
% a, b [mm]     : parallel faces of the quadrilateral
% c_1, c_2 [mm] : inclined segments of the quadrilateral
% h [mm]        : quadrilateral heigth
% r_1, r_2 [mm] : radii that define the curved faces of the quadrilateral
% L [mm]        : compressor axial length
%
%OUTPUT
% Vol [mm^3]    : basic volumes 
%
%NOTES
% see D Zanchi thesis for references and pictures (Capitolo 2)

    %% CALCULATIONS %%
    Area_trap = (a + b).*h*0.5;
    Area_circ_1 = (asin(0.5*c_1./r_1)).*r_1.^2 - 0.5*c_1.*(r_1.*cos(asin(0.5*c_1./r_1)));
    Area_circ_2 = (asin(0.5*c_2./r_2)).*r_2.^2 - 0.5*c_2.*(r_2.*cos(asin(0.5*c_2./r_2)));
    Area_quad = Area_trap - Area_circ_1 + Area_circ_2;
    Vol = Area_quad*L;
end
