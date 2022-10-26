function [mass_vane,moment_inertia_vane,theta_v,Y_F,DELTA,Sigma,LAMBDA,XI,Xg_vane,Yg_vane,b,r_tip]...
    = S2_Vane_Geom_Prop(b,r_tip,s,s_1D,l,rho_vane,offsetG,L,toll_t)
%this funcion calculates all the geometrical properties of a vane, such as characteristic
% angles, areas, total volume, mass, center of mass, moment of inertia.
%INPUT
% b                             [m]         : center of tip radius of curvature offset
% r_tip                         [m]         : tip radius
% s                             [m]         : vane thickness (2D geometry)
% s_1D                          [m]         : vane thickness (1D geometry)
% l                             [m]         : vane length
% rho_vane                      [m]         : vane density
% offsetG                       [%]         : offset vane mass center respect to mid plane
% L                             [m]         : axial length
% toll_t                        [-]         : teorical tollerance      
%
%OUTPUT
% mass_vane                     [kg]        : vane mass
% moment_inertia_vane           [kg*m^2]    : vane moment of inertia
% Y_F                           [m]         : y coordinates of center of curvature of tip
% Xg_vane, Yg_vane              [m]         : coordinates of center of mass
% theta_v, DELTA, Sigma,
% LAMBDA, XI                    [rad]       : characteristic angles of vane
%
%NOTES
% see D Zanchi thesis for references and pictures (Chapter 2)
% HISTORY: V10.1 DeFranco_Genoni_Gianoncelli: Area(1) modified to adjust areas evaluation in some peculiar cases 

    %% DECLARATIONS AND CALCULATIONS %%
    %In 1D geometry r_tip & b must be equal to 0 in order for the calculation to be correct
    if s==0 && r_tip ~= 0
        r_tip = 0;
        warning ('Make sure that r_tip = 0 for 1D geometry')
    end
    
    if s==0 && b ~= 0
        b = 0;
        warning('Make sure that b = 0 for 1D geometry')
    end

    % vectors preallocation of vane elementary zone
    Area = NaN(1,4);                                                % elemantary area
    Xg   = NaN(1,4);                                                % x-coord of center of mass
    Yg   = NaN(1,4);                                                % y-coord of center of mass

    % calculation of blade basic geometric elements - they are the same for every blade
    X_F     = 0.5*s - b;                                            % x-coord. of tip center of curvature
    
    if s == 0
       DELTA  = 0;                                                  % DEF angle (1D geometry)
       LAMBDA = pi/2;                                               % LFQ angle (1D geometry)
       XI     = pi/2;                                               % IFQ angle
    else
       DELTA  = asin(X_F/r_tip);                                    % DEF angle
       LAMBDA = acos((s-X_F)/r_tip);                                % LFQ angle
       XI     = acos((0.5*s-X_F)/r_tip);                            % IFQ angle
    end
    
    PI      = pi*0.5 - XI;                                          % GFI angle
    SIGMA   = XI - LAMBDA;                                          % ILF angle
    Sigma   = 0.5*(SIGMA + PI + DELTA);                             % LFE/2 angle (half of the tip sweep angle)
    theta_v = 0.5*pi - LAMBDA - Sigma;                              % 

    % definition of y-coord. of tip center of curvature - it depends on the x-coord position
    if X_F < 0 || X_F > s         
        Y_F = l-r_tip*cos(DELTA);
    else
        Y_F = l-r_tip;
    end

    % calculation of elementary areas necessary for center of mass evaluation
    Area(1) = s*Y_F;                                                %rectangle DQVZ
    Area(2) = X_F*r_tip*cos(DELTA)*0.5;                             %triangle EDF
    Area(3) = (s - X_F)*r_tip*sin(LAMBDA)*0.5;                      %triange LFQ
    Area(4) = Sigma*r_tip^2;                                        %circle sector EFL

    % calcuation of x-coord. elementary areas centers of mass
    Xg(1) = 0.5*s;                                                  % rectangle DQVZ
    Xg(2) = X_F/3;                                                  % triangle EDF
    Xg(3) = X_F + (2/3)*(s-X_F);                                    % triangle LQF

    if s == 0
        Xg(4) = 0;
    else
        Xg(4) = X_F + (2/3)*r_tip*sin(Sigma)*cos(Sigma+LAMBDA)/Sigma;   %x-coord. of circle sector EFL
    end
    
    % calcuation of x-coord. elementary areas centers of mass
    Yg(1) = 0.5*Y_F;                                                % rectangle DQVZ
    Yg(2) = Y_F + r_tip*cos(DELTA)/3;                               % rectangle EDF
    Yg(3) = Y_F + (1/3)*r_tip*sin(LAMBDA);                          % triangle LQF

    if s == 0
        Yg(4) = l;
    else
        Yg(4) = Y_F + (2/3)*r_tip*sin(Sigma)*sin(Sigma+LAMBDA)/Sigma;   %y-coord. of circle sector EFL
    end

    % calculation of vane center of mass coordinates via weighted average of the mass of each area
    
    if s == 0
        Xg_vane = 0;
        Yg_vane = Yg(4)/2 * (1 + offsetG/100);
    else
        Xg_vane = (Area*Xg')/sum(Area);
        Yg_vane = (Area*Yg')/sum(Area) * (1 + offsetG/100);         %the offset of vane G (if present) is considered
    end

    % calculation of vane mass and moment of inertia. Vane mass depends on
    % the type of simulation: if s and r_tip are zero, the vane model is 1D.
    if s < toll_t
        mass_vane = s_1D*l*L*rho_vane;                              %[kg]
    else
        mass_vane = sum(Area)*L*rho_vane;                           %[kg]
    end

    moment_inertia_vane = mass_vane*(l^2)/12;                       %[kg.m^2] computed as beam
    
    %% CHECKS %%
        % center of mass shall fall within the vane
        chck_1 = (Xg_vane < 0 || Xg_vane > s) || (Yg_vane > l || Yg_vane < 0);
        % checks on center of mass position
        chck_2 = Yg_vane > 0.5*l;
        chck_3 = b>0 && Xg_vane>0.5*s || b<0 && Xg_vane<0.5*s;

        if  chck_1 == 1
            warning('S2_Vane_Geom_Prop:barycenter','Center of mass does not fall within the vane')
            SX_Logfile ('w',{lastwarn});
        end

        if chck_2 ==1
            warning('S2_Vane_Geom_Prop:G_pos','Center of mass position is not compatible with vane geometry')
            SX_Logfile ('w',{lastwarn});
        end

        if chck_3 == 1
            warning('S2_Vane_Geom_Prop:G_pos','Center of mass position is not compatible with vane geometry')
            SX_Logfile ('w',{lastwarn});
        end
        
        if license('test','Symbolic_Toolbox')
            if contains(struct2array(ver), 'Symbolic Math Toolbox')
                syms limite;            
                chck_4 = double(limit(asin(X_F/limite), limite, r_tip));
                chck_5 = double(limit(acos((s-X_F)/limite), limite, r_tip));
                chck_6 = double(limit(acos((0.5*s-X_F)/limite), limite, r_tip));
                if abs(chck_4-DELTA) >toll_t || abs(chck_5-LAMBDA)>toll_t || abs(chck_6-XI)>toll_t
                    warning ('Incorrect geometry: check angles DELTA LAMBDA XI')
                end
                chck_7 = abs(Xg_vane - double(limit((limite*(l-r_tip)*Xg(1)+Area(2)*Xg(2)+Area(3)*Xg(3)+Area(4)*Xg(4))/(limite*(l-r_tip)+Area(2)+Area(3)+Area(4)),limite,s))) < toll_t;
                chck_8 = abs(Yg_vane - double(limit((limite*(l-r_tip)*Yg(1)+Area(2)*Yg(2)+Area(3)*Yg(3)+Area(4)*Yg(4))/(limite*(l-r_tip)+Area(2)+Area(3)+Area(4))*(1+offsetG/100),limite,s))) < toll_t;
                if chck_7 ~= 1 || chck_8 ~= 1
                    warning('Incorrect geometry: check center of mass Xg_vane Yg_vane')
                end
            else
                warning ('Download Symbolic Math Toolbox for geometry checks')
            end
        end

        clear chck_1 chck_2 chck_3 chck_4 chck_5 chck_6 chck_7 chck_8
    
    clear Area Xg Yg
end