function [rho_x,nu_x,mu_x,sigma_x,fOK] = S3_OilProperties(T_x,T_ref,rho_ref,nu_40,nu_100)
%This function calculates oil density, kinematic and dynamic viscosity and surface tension 
%at T_x
%  
%INPUT
% T_x [K]         : oil temperature
% T_ref [K]       : reference temperature for rho_ref
% rho_ref [kg/m3] : oil density @T_ref
% nu_40 [cSt]     : oil kinematic viscosity @40°C
% nu_100 [cSt]    : oil kinematic viscosity @100°C
%
%OUTPUT
% rho_x [kg/m3]   : oil density @T_x
% nu_x [cSt]      : oil kinematic viscosity @T_x
% mu_x [Pa*s]     : oil dynamic viscosity @T_x
% sigma_x [N/m]   : oil surface tension @T_x
% fOK [bool]      : OK flag (1 - no problem occourred, 0 - a problem have been spotted, SVEC will exit the simulation)
%
%
%EXAMPLES
%[rho_x,nu_x,mu_x,sigma_x] = S1_OilProperties(60,15,920,84,15.5);
%[rho_x,nu_x,mu_x,sigma_x] = S1_OilProperties([40 100],15,930,53,10.7);

  %% DEFINITIONS %%    
    fOK   = 1; 
    
  %% CONVERSIONS %%
    T_x   = T_x - 273.15;    % [K -> °C]
    T_ref = T_ref - 273.15;  % [K -> °C]
    
    %% DENSITY CALCULATIONS %%
    %Oil density @T_x is calculated by correlation proposed in "Estimation of
    %thermophisical properties of lubricating oils and their solutions with
    %refrigerants: an appraisal of existing methods" (Conde, 1996)
    A     = 0.6;
    rho_x = rho_ref-A*(T_x-T_ref); % [kg/m3]

    %% VISCOSITY CALCULATIONS %%
    %Kinematic
    nu_x = 10.^(log10(nu_40)+(T_x-40)/(100-40)*(log10(nu_100)-log10(nu_40))); %[cSt] 
    %Dynamic
    mu_x = nu_x/10^6.*rho_x; %[Pa*s]
    
    %% SURFACE TENSION CALCULATIONS %%
    %Tb (temperature at the normal boiling point of oil) and Tc (temperature at
    %the critical point of the oil) are evaluated using two equations simultaneously:
    %Nokay correlation (1959)         : Tc=15.2762*S^0.2985*Tb^0.62164
    %Kesler and Lee correlation (1976): Tc=189.8+450.6*S+(0.4244+0.1174*S)*Tb+(0.1441-1.0069*S)*10^5/Tb
    %The sistem made by these equations is solved within an interval of Tb value typical for 
    %Lubricant oil (400 K - 800 K) to avoid non physical solutions
    rho_H2Ol_60F = 999; % density of pure liquid water at 60°F (15.56°C) [kg/m3]
    rho_l_60F    = rho_ref-0.6*(15.56-T_ref); % oil density at 60°F (15.56°C)[kg/m3]
    S            = rho_l_60F/rho_H2Ol_60F;    % Specific gravity
    %search function
    fun = @(Tb) (189.8 +450.6 *S)+(0.4244+0.1174*S)*Tb +(0.1441-1.0069*S)*10^5/Tb -15.2762*S^0.2985*Tb^0.62164;
    %search procedure
    interval          = [400 800];
    [Tb,~,exitflag,~] = fzero(fun,interval); %Tb=[K]
    if exitflag~=1  
        warning( 'S1_OilProperties:NormalBoilingPoint','Oil normal boiling temperature not found')
        SX_Logfile ('w',{lastwarn});
        fOK = 0;
        Tb = NaN;
    end
    Tc = 15.2762*S^0.2985*Tb^0.62164; %[K]
    
    %pc (pressure at the critical point of the oil) is calculated by Kesler and Lee correlation (1976)
    pc =  (5.689   - 0.0566/S              ) ...
         -(0.43639 + 4.1216/S + 0.21343/S^2)/10^3 *Tb ...
         +(0.47579 + 1.182 /S + 0.15302/S^2)/10^6 *Tb^2 ...
         -(2.4505  +            9.9009 /S^2)/10^10*Tb^3 ; 
    pc = exp(pc); %[bar]
    %Surface tension @T_x is calculated by correlations proposed in "Surface tension
    %and the priciple of the corresponding states" (Brock and Bird, 1955)
    Q       = 0.1196*(1+(Tb/Tc*log(pc/1.01325))/(1-Tb/Tc))-0.279;
    sigma_x = pc^(2/3)*Tc^(1/3)*Q*(1-(T_x+273.15)/Tc).^(11/9); %[mN/m]
    sigma_x = sigma_x/1000; %[N/m]   
end

