function SVECpostProc(IO)
% This function manage the post processing of SVEC simulation results. The main features are:
% - Create, plot and eventually save figure
% - Plot and eventually save the report of the simulation
% All the saved files goes to the folder specified by IO struct
%
% INPUT
% IO  : Struct of basic input output parameters
%
% NOTE: * As it is defined now, figures size is approximately 16.5 x 17.5 cm (small change may occour)
%       * Figure size can be changed by modifing FigW and FigH

    %% PREAMBLE %%
    POSTtime = tic;                              % timer starts
    SX_Logfile('v',{'POSTprocessing started'});  % update logfile
    load (IO.fullloadname);                      % load simulation files
    
    %% UNPACK VARIABLES %%
    % unpack all variables with the coresponding name
    SX_v2struct(RESULTS);
    clear RESULTS

    %% PLOT CREATION %%
    if (IO.fDPLT || IO.fSPLT) && GEOMETRY.c
        % select figure size in centimeters
        FigW         = 17.5;                                            % figure width
        FigH         = 16.5;                                            % figure height
        OpenFig      = findobj(allchild(0),'flat', 'Type','figure');    % obtain list of opened figure
        [STDPlotSet] = S7_SetPlot (FigW,FigH,PROCESS.process,IO.fDPLT); % change plot properties for the current simulation
        
        % plot functions
        S7_PlotGeometry(PSI,Sigma,chck_tip,Npt);
        S7_PlotThermodynamics(theta,V_cell,V_comp,V_dis(1),p_cell,PROCESS.p_suc,p_out,T_g,T_mix,GAS.c_v,R_g,th_in,th_out,Gamma);
        S7_PlotNozzles(T_l,T_lmean,V_inj_cls,Dg,th_in,th_out,Gamma,NOZZLES)
        S7_PlotMechanics(F_m,F_v,F_t,T_m,T_v,T_t,F_cor,F_centr_i,F_centr_j,F_iz_i,conf,th_conf,Npt);
        S7_PlotPower (C1_pal_rot,C1_pal_P,C1_tip,C1_iz,C_izG,Npt);
        if FLAG.fLKG && FLAG.fLKG_plot
        S7_PlotLeakages(theta_vane,theta_SucOpen,theta_SucClose,theta_DisOpen,theta_DisClose,Flow_RS,Flow_VS,Flow_PE,Gamma);
        end
        if FLAG.fSTR
            S7_PlotStress(posMaxMf,N_posMaxMf,T_posMaxMf,Mf_posMaxMf,sVB,sBI,NUMERIC.stp,pos_contct);
        end
        
        % restore default matlab properties for plots
        reset (groot);
        propname = fieldnames(STDPlotSet);
        for iprop = 1:length(propname)
            set (groot, propname{iprop}, STDPlotSet.(propname{iprop}));
        end
    end

    %% PLOT SAVING %%
    if IO.fSPLT && GEOMETRY.c
        FigDir  = [IO.fullsavename ' Fig'];                              % name of figure directory
        mkdir     (FigDir);                                              % create figure directory
        FigList = findobj(allchild(0),'flat','Type', 'figure');          % obtain figure list
        for iFig = 1:length(FigList) - length(OpenFig)                   % for cycle to save each new figure
            FigHandle = FigList(iFig);                                   % fig index
            FigName   = get(FigHandle, 'Name');                          % get assigned name to save the figure
            set(FigHandle, 'CreateFcn', 'set(gcbo,''Visible'',''on'')'); % set figure as 'visible' when they will be opened again
            saveas(FigHandle, fullfile( FigDir,FigName), IO.ext);        % save figure in the specified flder, with the specified extention
        end
        InvFig = findobj(allchild(0),'flat','Type', 'figure','Visible','off');  % find invisible figure
        close (InvFig)                                                          % close invisible figure
    end

    %% REPORT %%
    if IO.fDRPT && GEOMETRY.c
        S7_PrintIO(1,IO,FLAG,NUMERIC,PROCESS,GEOMETRY,GAS,LIQUID,LEAK,VANE,SHAFT,BUSHING,NOZZLES,th_flp,th_Vmax,V_inj_nzl,rho_o,nu_o,sigma_o,Dg,fg,epsm,cls_indx,p_out,p_geom,T_g,T_lmean,T_mix,n_pltr,Q_th,L_sp_i,L_sp_m,E_sp_m,eta_g,eta_gT,eta_gl,eta_vol,eta_mecc,m_gas_Id,m_gas,q_gas_Id,q_gas,m_liq_Id,m_liq,Pind_PV,C_shaft_mean,P_mecc,M,Q_g_RS_suction,Q_g_PE_suction,Q_g_VS_suction,Q_g_PE_closed_cell,Q_g_VS_closed_cell,Fleak_g,Fleak_g_PE_suc,Fleak_g_VS_suc,Fleak_g_PE_cell_closed,Fleak_g_VS_cell_closed,Q_g_RS_dis,Q_g_PE_dis,Q_g_VS_dis)
        if FLAG.fSTR
            S7_PrintStress(1,IO,STRESS,STRESSsh,posMaxMf,Npt,MaxMf,Fs_traz,Fs_FA,sigmaI_max,aGi,aGcentr_i,aGcentr_j,aGcor_j,p_1,p_2,VN_i,VM_i,conf,FB,kB,pB,vB,pvB,sigmaeqD,rSFsD,rSFfD,xscrD,thcrD,sigmaeqDsh,rSFsDsh,rSFfDsh,xscrDsh,thcrDsh);
        end
    end
    if IO.fSRPT && GEOMETRY.c
        fileID = fopen ([IO.fullsavename,' report.txt'],'wt');
        S7_PrintIO(fileID,IO,FLAG,NUMERIC,PROCESS,GEOMETRY,GAS,LIQUID,LEAK,VANE,SHAFT,BUSHING,NOZZLES,th_flp,th_Vmax,V_inj_nzl,rho_o,nu_o,sigma_o,Dg,fg,epsm,cls_indx,p_out,p_geom,T_g,T_lmean,T_mix,n_pltr,Q_th,L_sp_i,L_sp_m,E_sp_m,eta_g,eta_gT,eta_gl,eta_vol,eta_mecc,m_gas_Id,m_gas,q_gas_Id,q_gas,m_liq_Id,m_liq,Pind_PV,C_shaft_mean,P_mecc,M,Q_g_RS_suction,Q_g_PE_suction,Q_g_VS_suction,Q_g_PE_closed_cell,Q_g_VS_closed_cell,Fleak_g,Fleak_g_PE_suc,Fleak_g_VS_suc,Fleak_g_PE_cell_closed,Fleak_g_VS_cell_closed,Q_g_RS_dis,Q_g_PE_dis,Q_g_VS_dis)
        if FLAG.fSTR
            S7_PrintStress(fileID,IO,STRESS,STRESSsh,posMaxMf,Npt,MaxMf,Fs_traz,Fs_FA,sigmaI_max,aGi,aGcentr_i,aGcentr_j,aGcor_j,p_1,p_2,VN_i,VM_i,conf,FB,kB,pB,vB,pvB,sigmaeqD,rSFsD,rSFfD,xscrD,thcrD,sigmaeqDsh,rSFsDsh,rSFfDsh,xscrDsh,thcrDsh);
        end
        fclose(fileID);
    end
    
    SX_Logfile('v',{'POSTprocessing completed. Elapsed time: %.2f s ',toc(POSTtime)});
end