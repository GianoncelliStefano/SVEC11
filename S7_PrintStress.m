function S7_PrintStress(fileID,IO,STRESS,STRESSsh,posMaxMf,Npt,MaxMf,Fs_traz,Fs_FA,sigmaI_max,aGi,aGcentr_i,aGcentr_j,aGcor_j,p_1,p_2,VN_i,VM_i,conf,FB,kB,pB,vB,pvB,sigmaeqD,rSFsD,rSFfD,xscrD,thcrD,sigmaeqDsh,rSFsDsh,rSFfDsh,xscrDsh,thcrDsh)
% This function creates and plot a report with the user input and main
% simulation reults concerning the stress analysis
%  INPUT
% fileID [-]        : File identifier of the opened file
% IO [-]            : structure with basic input-output parameters 
% STRESS [-]        : Structure of stress parameters
% STRESSsh [-]      : Structure of shaft stress parameters
% posMaxMf [rad]    : discretized angular position where maxMf has been found
% Npt [-]           : number of grid discretization points in [0:2pi]
% MaxMf [Nm]        : absolute maximum of Mf
% Fs_traz [-]       : static tensile safety factor
% Fs_FA [-]         : fatigue safety factor
% sigmaI_max [pa]   : greatest I main stress obtained via Mohr circle
% aGi [m/s2]        : axial component of vane center of mass inertial acceleration (analitical)
% aGcentr_i [m/s2]  : axial component of vane center of mass centripetal acceleration (axial to vane)
% aGcentr_j [m/s2]  : perpendicular component of vane center of mass centripetal acceleration (perpendicular to vane)
% aGcor_j [m/s2]    : perpendicular component of vane center of mass Coriolis acceleration (perpendicular to vane)
% p_1 [pa]          : low pressure acting on vane in [0:2pi]
% p_2 [pa]          : high pressure acting on vane in [0:2pi]
% VN_i [m]          : vane side surfaces inside the rotor slot - suction side
% VM_i [m]          : vane side surfaces inside the rotor slot - pressure side
% conf [-]          : array of vane configurations
% FB  [N]                : maximum bushing load
% kB  [deg]              : bushing load direction
% pB  [N/mm^2]           : bushing specific load
% vB  [m/s]              : sliding velocity on bushings
% pvB [N/mm^2 • m/s]     : bushing Pressure-Volume Value
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

    %% DEFINITIONS %%
    % Titles
    title0 = '****************** STRESS ANALYSIS  ******************';
    title1 = 'STRESS PARAMETERS';
    title2 = 'SHAFT STRESS ANALYSIS';
    title3 = 'RESULTS: TENSILE STRESS';
    title4 = 'STATIC CHECK PARAMETERS FOR FEM ANALYSIS';
    title5 = 'RESULTS: SHAFT + ROTOR STRESS';
    title6 = 'RESULTS: SHAFT STRESS';
    title7 = 'RESULTS: BUSHING PERFORMANCES';

    % Divider
    sep1   = '******************************************************';
    sep2   = '------------------------------------------------------';

    % Output formats
    f0 = '%s \n';
    f2 = '%-32s %-10s %-15.1f \n';
    f4 = '%-32s %-10s %-15.2f \n';
    f7 = '%-32s %-10s %-15.0f \n';
    f8 = '%-32s %-10s %-15.3f \n';

    % Unità di misura
    adim = '-';
    gradi = '°';
    mm = 'mm';
    bar = 'bar';
    acc = 'm/s2';

    %% PLOT %%
    % Plot report title
    fprintf(fileID,'\n%s\n',sep1);
    fprintf(fileID,f0,title0);
    fprintf(fileID,'%s',sep1);

    % Plot cast iron mechanical properties
    fprintf(fileID,'\n%s \n',sep2);
    fprintf(fileID,f0,title1);
    fprintf(fileID,'%s \n',sep2);
    fprintf(fileID,f7,'Ultimate tensile stress','MPa',STRESS.sigmaR_traz*1e-6);
    fprintf(fileID,f7,'Ultimate compressive stress','MPa',STRESS.sigmaR_comp*1e-6);
    fprintf(fileID,f4,'Fatigue stress ratio',adim,STRESS.ratio_FA);
    fprintf(fileID,f4,'Fatigue size factor',adim,STRESS.b2);
    fprintf(fileID,f4,'Fatigue surface finish factor',adim,STRESS.b3);
    
    % =====================================================================
    % =====================================================================
    fprintf(fileID,'\n%s \n',sep2);
    fprintf(fileID,f0,title2);
    fprintf(fileID,'%s \n',sep2);
    fprintf(fileID,f7,'Ultimate tensile stress','MPa',STRESSsh.Rmt_s*1e-6);
    fprintf(fileID,f7,'Ultimate compressive stress','MPa',STRESSsh.Rmc_s*1e-6);
    fprintf(fileID,f7,'Yield stress','MPa',STRESSsh.Rp02_s*1e-6);
    fprintf(fileID,f7,'Fatigue stress','MPa',STRESSsh.sigma_fas*1e-6);
    fprintf(fileID,f4,'Fatigue size factor',adim,STRESSsh.b2_s);
    fprintf(fileID,f4,'Fatigue surface finish factor',adim,STRESSsh.b3_s);
    % =====================================================================
    % ===================================================================== 

    % Plot tensile stress analysis
    fprintf(fileID,'\n%s \n',sep2);
    fprintf(fileID,f0,title3);
    fprintf(fileID,'%s \n',sep2);
    fprintf(fileID,f2,'Vane angular position',gradi,posMaxMf*360/Npt);
    fprintf(fileID,f2,'Maximum bending action','N m',MaxMf);
    fprintf(fileID,f2,'I principal stress','MPa',sigmaI_max*1e-6);
    fprintf(fileID,f2,'Tensile safety coefficient',adim,Fs_traz);
    fprintf(fileID,f2,'Fatigue safety coefficient',adim,Fs_FA);

    % Plot static check parmeter for FEM analysis
    fprintf(fileID,'\n%s \n',sep2);
    fprintf(fileID,f0,title4);
    fprintf(fileID,'%s \n',sep2);
    fprintf(fileID,f2,'Pressure side 1',bar,p_1(posMaxMf)*1e-5);
    fprintf(fileID,f2,'Pressure side 2',bar,p_2(posMaxMf)*1e-5);
    fprintf(fileID,f2,'Side 1 concealed surface (VN)',mm,VN_i(posMaxMf)*1e+3);
    fprintf(fileID,f2,'Side 2 concealed surface (VM)',mm,VM_i(posMaxMf)*1e+3);
    fprintf(fileID,f2,'Inertial acceleration',acc,aGi(posMaxMf));
    fprintf(fileID,f2,'Radial centr. acceleration',acc,aGcentr_i(posMaxMf));
    fprintf(fileID,f2,'Tangential centr. acceleration',acc,aGcentr_j);
    fprintf(fileID,f2,'Coriolis acceleration',acc,aGcor_j(posMaxMf));
    fprintf(fileID,f2,'Vane configuration',adim,conf(posMaxMf));
    
    % =====================================================================
    % ===================================================================== 
    % Plot shaft+rotor stress analysis
    fprintf(fileID,'\n%s \n',sep2);
    fprintf(fileID,f0,title5);
    fprintf(fileID,'%s \n',sep2);
    fprintf(fileID,f4,'Vane angular position',gradi,thcrD*180/pi);
    fprintf(fileID,f2,'Position on the shaft+rotor',mm,xscrD*1e3);
    fprintf(fileID,f4,'Maximum shaft+rotor stress','MPa',sigmaeqD);
    fprintf(fileID,f8,'Ductile static safety factor',adim,rSFsD);
    fprintf(fileID,f8,'Ductile fatigue safety factor',adim,rSFfD);
    
    % Plot shaft stress analysis
    fprintf(fileID,'\n%s \n',sep2);
    fprintf(fileID,f0,title6);
    fprintf(fileID,'%s \n',sep2);
    fprintf(fileID,f4,'Vane angular position',gradi,thcrDsh*180/pi);
    fprintf(fileID,f2,'Position on the shaft+rotor',mm,xscrDsh*1e3);
    fprintf(fileID,f4,'Maximum shaft stress','MPa',sigmaeqDsh);
    fprintf(fileID,f8,'Ductile static safety factor',adim,rSFsDsh);
    fprintf(fileID,f8,'Ductile fatigue safety factor',adim,rSFfDsh);
    
    % Plot bushing performance analysis
    fprintf(fileID,'\n%s \n',sep2);
    fprintf(fileID,f0,title7);
    fprintf(fileID,'%s \n',sep2);
    fprintf(fileID,f8,'Maximum bushing load','kN',FB/1e3);
    fprintf(fileID,f4,'Bushing load direction','deg',mean(kB));
    fprintf(fileID,f4,'Maximum specific bushing load','N/mm2',pB);
    fprintf(fileID,f4,'Sliding velocity on bushings','m/s',vB);
    fprintf(fileID,f4,'Bushing Pressure-Volume value','N/mm2*m/s',pvB);
    fprintf(fileID,'\n%s \n',sep2);
    % =====================================================================
    % ===================================================================== 
      
    fprintf(fileID, 'Elapsed time: %.2f s \n',toc(IO.CALLtime));
end