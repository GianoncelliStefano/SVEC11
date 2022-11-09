 function [fOK] = SVECmodelMain(IO)
% This function is divided in 6 section:
% Section S1: take input from structure memorised in CALL_SVECmodelMain and perform input checking and conversion
% Section S2: solve geometric model
% Section S3: solve thermodynamic model
% Section S4: solve dynamic model
% Section S5: solve efficiency and power
% Section S6: solve internal actions model

    %% PREAMBLE %%
    MAINtime = tic;
    SX_Logfile('v',{'SVECmodelMAIN started'}); % update logfile   
    load  (IO.fullloadname);                   % load input structures
    
    %% S1 - INPUT CHECK & UNPACKING %%
    [fOK] = S1_InputCheck(IO,FLAG, NUMERIC, PROCESS, GEOMETRY, GAS, LEAK, VANE, NOZZLES);
    
   
    % Unpack variables
    SX_v2struct(FLAG);
    SX_v2struct(NUMERIC);
    SX_v2struct(PROCESS);
    SX_v2struct(GEOMETRY);
    SX_v2struct(GAS);
    SX_v2struct(LIQUID);
    SX_v2struct(VANE);
    SX_v2struct(STRESS);
    SX_v2struct(LEAK);
    SX_v2struct(STRESSsh);
    SX_v2struct(SHAFT);
    SX_v2struct(BUSHING);
    SX_v2struct(INTAKE);

    % ================================================================
    % ================================================================
    % Addition of new structuresas input for SVECpreProc
    % ================================================================
    % ================================================================
    
    %% S2 - GEOMETRY & KINEMATICS %%
    % NB: For geometric computations, the reference is the TRAILING VANE

    % Grid Discretization
    if fOK, [Npt,Npt_cell,Gamma,theta, pos_tang, theta_vane] = S2_Discretization(Npt_i,n_van,c);
    end
    
    % ===============================================================
    % ===============================================================
    % Discretization of the shaft for internal action analysis
    if fOK, [xsv] = S2_ShaftDiscretization(bB,d1,L,shm);
    end
    % ===============================================================
    % ===============================================================

    if c == 1    % Circular stator
        % calculation of geometric properties of vane
        if fOK, [mass_vane,moment_inertia_vane,theta_v,Y_F,DELTA,Sigma,LAMBDA,XI,Xg_vane,Yg_vane]...
                = S2_Vane_Geom_Prop(b, r_tip, s,s_1D, l , rho_vane, offsetG, L, toll_t);
        end

        % Calculation of geometric properties of machine
        if fOK, [Rth,PSI,BI,HI,UPSILON,UPSILON_tip_1,UPSILON_tip_2,ME,NL,EP_arch,PL_arch,BM,BN,AX_i,AY_i,AP_i,AP_j,EI,IL,AI,AI_1D,AR_i,AR_j,AS_i,AS_j,AV_i,AM_i,R_v,VN_i,VM_i,sgn]...
                = S2_Geometry_C(D,d,theta_v,theta_vane,th_tilt,Sigma,r_tip,XI ,b,s,DELTA,LAMBDA,toll_t,Y_F);
        end

        % Calculation of kinematics quantities
        if fOK, [AG,b_gx,xG,rv,psi,alpha,vAGi,aGi,aGcentr_i,aGcentr_j,aGcor_j,aIVcentr_i,omega,DomegaG]...
                = S2_Kinematics_C(D,d,b,Y_F,Xg_vane,Yg_vane,R_v,theta_vane,rpm,r_tip,s,l,th_tilt,sgn,AV_i,AI,Npt,toll_t,toll_d,fDBG);
        end

        % Calculation of volumes to be excluded
        if fOK, [Vol_Excluded] = S2_Volume_Excluded_C( BI, BN, BM, ME, NL, EI, IL, theta_vane, L, s, r_tip, d, AI_1D, AI, AP_i, AP_j, D, toll_t,fDBG);
        end

        % Calculation of volume between two vanes during a complete rotor revolution [m^3]
        if fOK, [V_cell] = S2_Volume_C(D,d,L,theta,Gamma,sgn,Npt_cell,pos_tang,th_tilt,toll_d,fDBG,Vol_Excluded,s);
        end

    elseif c == 2     % Eliptical stator 
        % Geometric and kinematics quantities calculations
        if fOK, [mass_vane,Rth,xG,w,beta,psi,vGr,aGr_iz,aGr_centr,aGt_cor,omega,sgn] = S2_Kinematics_E(d,e,s,L,rho_vane,theta_vane,rpm,l,toll_d,Npt,fDBG);
        end

        % Calculation of volume between two vanes during a complete rotor revolution [m^3]
        if fOK, [V_cell] = S2_Volume_E(d,e,L,theta,Gamma,Npt_cell,pos_tang);
        end
    end

    % Ports angles and positions needs to be corrected
    if fOK, [pos_SucOpen,pos_SucClose,pos_DisOpen,pos_DisClose,pos_End,th_in,th_out,th_flp,th_Vmax,fOK] = S2_Angles (c,n_van,process,StdProcess,theta,th_SucOpen,th_SucClose,th_DisOpen,th_DisClose,TgAngle,V_cell,pos_tang,sgn,psi,Gamma,NOZZLES);
    end

    % cell volume during suction, compression and delivery processes [m^3]
    if fOK
        [V_suc,V_dis,V_comp] = S2_VolumeSplit (V_cell,pos_SucClose,pos_DisOpen,Npt_cell);
    end

    %% SINLET- INLET DUCT+PORT MODEL %%
    % Inlet process model. Set fSDP to 0 for air-end only computation
    if fSDP
     if fOK, [p_suc,T_suc] = S2BIS_inlet(p_suc , T_suc , INport_Amax, INport_Amin, V_comp(1), MM_g , n_van,rpm,c,toll_d,c_v,coeff_invalve,pipe,cpitch,ct,lenght,D_up,D_do,roughness,mu_g,coeff_infilter);
     end
    end
    
    %% S3 - THERMODYNAMICS %%
    % Computation of pressure and temperature during closed-cell process
    if fOK, [Q_g_RS_dis,Q_g_PE_dis,Q_g_VS_dis,Q_g_RS_suction,Q_g_PE_suction,Q_g_VS_suction,Q_g_PE_closed_cell,Q_g_VS_closed_cell,Fleak_g,Fleak_g_PE_suc,Fleak_g_VS_suc,Fleak_g_PE_cell_closed,Fleak_g_VS_cell_closed,m_g,q_gas_Id,q_gas,T_l,T_g,T_mix,T_lmean,p_comp,rho_o,nu_o,sigma_o,m_inj_nzl,V_inj_nzl,V_inj_cls,epsm,W_cell,n_pltr,Dg,fg,cls_indx,M,m_gas_Id,m_gas,m_liq_Id,m_liq,Q_th,R_g,fOK,theta_vane,theta_SucOpen,theta_SucClose,theta_DisOpen,theta_DisClose,Flow_RS,Flow_VS,Flow_PE] ...
            = S3_Thermodynamics(process,D,d,d_hub,L,BI,s,RSclr,VSclr,PEclr,TgAngle,NOZZLES,rho_ref,T_ref,nu_40,nu_100,c_l,c_v,k_g,k_l,MM_g,T_suc,T_0,p_0,p_suc,p_del,mu_g,V_cell,V_comp,pos_SucOpen,pos_SucClose,pos_DisOpen,pos_DisClose,pos_End,Npt_cell,theta,theta_vane,omega,toll_d,Gamma,c,n_van,rpm,model_g,model_full,name4prop,molcomp_g,propSource,MUmodel_RS,MUmodel_VS,MUmodel_PE,RSmodel,VSmodel,PEmodel,IterMax,UndRlx,fLKG,fLKG_in);
    end
        
    % Different pressure vectors are computed
    if fOK, [p_out,p_cell,p_1,p_2,p_geom] = S3_Pressure(theta,theta_vane,p_comp,pos_SucClose,pos_DisOpen,Npt_cell,p_del);
    end

%       %% SOUTLET - OUTLET DUCT+PORT MODEL %%
%   Outlet process model. Set fSDP to 0 for air-end only computation
%   if fSDP
%     if fOK, [p_u,T_u] = Soutlet(p_del , T_mix(end) , OUTport_Amax, OUTport_Amin,m_gas, MM_g, c_v, p_geom);
%     end
%   end

    %% S4 - MECHANICS %%
    if c == 1
        % Calculation of forces and acting on a single vane [NOTE: fix the checks]
        if fOK, [F_t,T_t,F_m,T_m,F_v,T_v,F_cor,F_iz_i,F_centr_i,F_centr_j,Fw_i,Fw_j,Fp_1,Fp_2,Fp_1_tip,Fp_2_tip,Fp_cava,Fp_1_cava,Fp_2_cava,C_izG,AM_i_new,Ay_i,Ax_i,Mf_ax,conf,th_conf,s1,s2,chck_tip]...
                = S4_VaneDynamics_C(f_c,f_t,L,s,th_tilt,theta_vane,mass_vane,moment_inertia_vane,AG,b_gx,UPSILON,UPSILON_tip_1,UPSILON_tip_2,ME,NL,EP_arch,PL_arch,AX_i,AY_i,AP_j,AP_i,AR_j,AR_i,AS_j,AS_i,AV_i,AM_i,VN_i,VM_i,vAGi,aGcentr_i,aGcentr_j,aGi,aGcor_j,DomegaG,sgn,p_1,p_2,Npt,rot_dir,zeta,toll_t,fDBG);
        end

        % Calculation of torque acting on single vane and on shaft [NOTE: fix the checks]
        if fOK, [C1_pal_rot,C1_pal_P,C1_tip,C1_iz] = ...
            S4_VaneTorque_C(AG,b_gx,s,rv,F_t,T_t,F_m,T_m,F_v,T_v,F_cor,F_centr_i,F_centr_j,Fw_i,Fw_j,F_iz_i,Fp_1,Fp_2,Fp_1_tip,Fp_2_tip,Fp_cava,Fp_1_cava,Fp_2_cava,C_izG,UPSILON,UPSILON_tip_1,UPSILON_tip_2,AX_i,AY_i,AP_i,AP_j,AR_i,AR_j,AS_i,AS_j,AV_i,AM_i_new,Ay_i,Ax_i,sgn,toll_t,s1,s2);
        end

        % Shaft reaction forces (necessary for friction in bushing) [NOTE: add rotor weight force]
        if fOK, [R1_x,R1_y] = S4_VaneRotorForces_C(d,L,n_van,theta_vane,alpha,BM,BN,F_m,F_v,T_m,T_v,Fp_cava,Fp_1_cava,Fp_2_cava,p_cell,c,Npt,Npt_cell,sgn,th_tilt);
        end

    elseif c==2
        % Calculation of forces acting on a single vane
        if fOK, [F_t,T_t,F_m,T_m,F_v,T_v,F_p,F_centr_i,F_iz_i,F_cor,conf,chck_tip] = S4_VaneDynamics_E(f_c,f_t,L,l,w,Rth,d,xG,beta,mass_vane,vGr,aGr_centr,aGr_iz,aGt_cor,p_1,p_2,Npt,toll_t,fDBG);
        end

        % Calculation of torque acting on single vane and on shaft
        if fOK, [C1_pal_rot,C1_pal_P,C1_tip,C1_iz] = S4_VaneTorque_E(Rth,d,xG,l,w,beta,F_t,F_m,F_v,T_t,F_p,F_cor,toll_t);
        end

        % Shaft reaction forces (necessary for friction in bushing)
        if fOK, [R1_x,R1_y] = S4_VaneRotorForces_E(d,L,n_van,theta_vane,F_m,F_v,T_m,T_v,p_cell,c,Npt,Npt_cell);
        end
    end
    
    % =====================================================================
    % =====================================================================
    % Rotor+shaft weight force calculation
    if fOK, [Fg,qG] = S4_RotWeight(L,d,d_hub,l,s,d1,bB,n_van,rho_s,xsv);
    end
        
    % Computation of overall torque on shaft
    if fOK, [R_x,R_xg,R_y,R_yg,R_rot,C_rot,C_bronz,C_shaft,C_press,C_shaft_mean] = S4_ShaftTorque(C1_pal_rot,C1_pal_P,R1_x,R1_y,n_van,Npt_cell,d_hub,f_b,Fg,zeta,rot_dir,c);
    end
    
    % =====================================================================
    % =====================================================================

    %% S5 - POWER & EFFICIENCY %%
    % Computes indicated power e work performed by thermodynamic model
    if fOK, [Wtot,Wis_g,Wis_gl,Wis_T,Pind_PV] = S5_PvPower(W_cell,p_suc,p_out,V_suc,V_dis,rpm,n_van,c,R_g,c_v,c_l,T_suc,m_g,m_inj_nzl,V_inj_nzl);
    end

    % Computes indicated power performed by mechanic model
    if fOK, [P_mecc,a] = S5_MechPower(C_shaft,C_press,omega,Pind_PV,toll_d);
    end

    % Computes thermodynamic and mechanic efficiencies
    if fOK, [eta_g,eta_gl,eta_gT,eta_vol,eta_mecc,L_sp_m,L_sp_i,E_sp_m] = S5_Efficiency(Wtot,Wis_g,Wis_gl,Wis_T,Pind_PV,P_mecc,m_gas_Id,m_gas,q_gas,m_inj_nzl,toll_d,process);
    end

    %% S6 - INTERNAL ACTIONS %%
    if fSTR && c==1
        % 1D geometry segments used for internal action computation
        if fOK, [VB,VI] = S6_ApplicationPointForces(Y_F,HI, AV_i, d, th_tilt);
        end
        
        % Discretization of vane for internal action analysis
        if fOK, [sVB,sBI,pos_contct] = S6_SaxisVane(BI,VB,stp,Npt);
        end
        
        % Internal action calculation
        if fOK, [MaxMf,N_sMaxMf,T_sMaxMf,posMaxMf,N_posMaxMf,T_posMaxMf,Mf_posMaxMf] = S6_InternalActions(L,s,mass_vane,UPSILON,UPSILON_tip_1,UPSILON_tip_2,VI,sVB,sBI,aGi,aIVcentr_i,F_v,T_v,F_m,T_m,F_t,T_t,Fp_cava,Fp_1_tip,Fp_2_tip,F_cor,F_centr_j,Fw_i,Fw_j,C_izG,Mf_ax,Npt,p_2,p_1,pos_contct,sgn,toll_d,toll_t);
        end
        
        % Safety coefficent computation for most strained point
        if fOK, [sigmaI_max,Fs_traz,Fs_FA] = S6_StateOfStress(MaxMf,N_sMaxMf,T_sMaxMf,s_1D,L,sigmaR_traz,sigmaR_comp,ratio_FA,b2,b3);
        end
   % ====================================================================
   % ====================================================================
        % Bushing mechanical performance calculation
        if fOK, [FB,kB,pB,vB,pvB] = S6_BushingForces(bB,R_rot,R_xg,R_yg,rpm,d_hub,zeta,rot_dir);
        end
   
        % Evaluation of shaft internal action
        if fOK, [thv_res,Tx,Ty,Mx,My,Mz,pck1,pck2,pck3,pck4,pck5] = S6_RotAct(xsv,L,d1,bB,R_x,R_xg,R_y,R_yg,qG,C_bronz,C_rot,theta,Npt_cell,shm,zeta,rot_dir,toll_t,n_van,pos_SucOpen,pos_SucClose,pos_DisOpen,pos_DisClose);
        end
        
        % Rotor+shaft static and fatigue stress analysis
        if fOK, [sigmaeqD,rSFsD,rSFfD,xscrD,thcrD,sigmaeqDsh,rSFsDsh,rSFfDsh,xscrDsh,thcrDsh] = S6_RotStress(thv_res,xsv,pck1,pck2,pck3,pck4,pck5,d_hub,d,l,Tx,Ty,Mx,My,Mz,Rmt_s,Rp02_s,sigma_fas,b2_s,b3_s);
        end
   % ====================================================================
   % ==================================================================== 
     
    end
    
    %% SAVE %%
    % SVEC pack and save in RESULTS structure only the data NOT contained in the input structure.
    % If add variables in this workspace and you need them to be post-processed, just addd their name here.         
    % Results are saved in two structure: 
    %  - RESULTS_A are data that are always generated (geometry, thermoduynamic, dynamic) 
    %  - RESULTS_B is for stress analysys, that may be not generated 
    
    if fOK && c == 1
        % Creation of result structures and saving
        RESULTS_A = SX_v2struct(Npt,Sigma,Gamma,PSI,theta,V_cell,V_comp,V_dis,p_cell,p_out,T_g,T_l,T_mix,V_inj_cls,Dg,R_g,th_in,th_out,th_flp,th_Vmax,...
                         F_m,F_v,F_t,T_t,T_m,T_v,F_cor,F_centr_i,F_centr_j,F_iz_i,conf,th_conf,chck_tip,C1_pal_rot,C1_pal_P,C1_tip,C1_iz,C_izG,Pind_PV,C_shaft_mean,P_mecc,...
                         rho_o,nu_o,sigma_o,fg,epsm,cls_indx,p_geom,T_lmean,V_inj_nzl,n_pltr,Q_th,L_sp_i,L_sp_m,E_sp_m,eta_g,eta_gT,eta_gl,eta_vol,eta_mecc,...
                         m_gas_Id,m_gas,q_gas_Id,q_gas,m_liq_Id,m_liq,M,aGi,aGcentr_i,aGcentr_j,aGcor_j,p_1,p_2,VN_i,VM_i,Fleak_g,Fleak_g_PE_suc,Fleak_g_VS_suc,Fleak_g_PE_cell_closed,Fleak_g_VS_cell_closed,s,...
                         Q_g_RS_dis,Q_g_PE_dis,Q_g_VS_dis,Q_g_RS_suction,Q_g_PE_suction,Q_g_VS_suction,Q_g_PE_closed_cell,Q_g_VS_closed_cell,theta_vane,theta_SucOpen,theta_SucClose,theta_DisOpen,theta_DisClose,Flow_RS,Flow_VS,Flow_PE);
        if fSTR
            RESULTS_B = SX_v2struct(posMaxMf,N_posMaxMf,T_posMaxMf,Mf_posMaxMf,sVB,sBI,pos_contct,MaxMf,Fs_traz,Fs_FA,sigmaI_max,...
                             FB,kB,pB,vB,pvB,sigmaeqD,rSFsD,rSFfD,xscrD,thcrD,sigmaeqDsh,rSFsDsh,rSFfDsh,xscrDsh,thcrDsh);
            RESULTS   = SX_Catstruct(RESULTS_A,RESULTS_B);
        else
            RESULTS   = RESULTS_A;
        end
              
        save (IO.fullsavename,'FLAG','NUMERIC','PROCESS','GEOMETRY','GAS','LIQUID','LEAK','VANE','SHAFT','BUSHING','NOZZLES','STRESS','STRESSsh','RESULTS');
        SX_Logfile('v',{'SVECmodelMAIN completed. Elapsed time: %.2f s ',toc(MAINtime)})
    else
        fprintf('\nSVEC_Main: an error occoured. SVEC could not complete the simulation. Read the previuos warning messages or the log file to try to fix the problem\n');
        SX_Logfile('a',{'SVECmodelMAIN aborted. Check the log file to try to fix the problem'})
    end

end
