function  [dOl_dp,dY_dp,dFc_dp,ddeltat_dp,der_HCcH,der_HCcHF,B12out,Hc,B] = derOl_Y_F_deltat_dpi(A,Zdes,angleact,Ccc,der_Cc_dp,ind_cont,ind_slip,Fc_mc,Fc_mc3,P_mc,Pfree_c3,contAngle,psi_first,lp,dPFreec_dp,dPcq_1_dp)
    if contAngle==1
        N_ini = gradient_simu_N_ini_contact(ind_cont,Pfree_c3,Fc_mc3,angleact,psi_first,contAngle);
        Y = gradient_simu_Y_contact(ind_cont,Rini,Fc_mc3,dPFreec_dp,lp,contAngle,[]);
        dOl_dY = (E * A_1B + N_ini) \ (EinvTheta + Y);
    else
        Y = gradient_simu_Y_contact(ind_cont,Rini,Fc_mc3,dPFreec_dp,lp,contAngle,dPcq_1_dp);
        dOl_dY = (E * A_1B) \ (EinvTheta + Y);
    end
    dF_ddelta_dp = A_1B * dOl_dY + inv_A_theta;
    dFc_dp = dF_ddelta_dp(1:3*length(ind_cont),:);
    ddeltat_dp = dF_ddelta_dp(3*length(ind_cont)+1:end,:);

    dOl_dp = dOl_dY(1:3,:);
    dY_dp = dOl_dY(4:6,:); 
end