function  [dOl_dp,dY_dp,dFc_dp,ddeltat_dp,der_HCcH,inv_A,A_1B,EA_1B,E] = derOl_Y_F_deltat_dp(A,Zdes,angleact,Ccc,der_Cc_dp,ind_cont,ind_slip,Fc_mc,Fc_mc3,P_mc,Pfree_c3,contAngle,psi_first,lp,dPFreec_dp,B_alg,dPcq_1_dp)
    theta = angleact(1);
    phi = angleact(2);
    psi = angleact(3);
    R = Rot(theta,phi,psi);
    Rini = Rot(theta,phi,psi_first);
    %B = gradient_simu_B_contact(angleact,Pfree_c3,Fc_mc3,ind_cont,Ccc,ind_slip,psi_first,contAngle);
    E = gradient_simu_E_contact(ind_cont,ind_slip,Fc_mc,P_mc,Zdes);
    [Theta_mat,der_HCcH] = gradient_simu_Theta_contact(ind_cont,der_Cc_dp,ind_slip,Fc_mc3,R,Rini,lp,dPFreec_dp,contAngle,dPcq_1_dp);
    inv_A = inv(A);
    A_1B = inv_A * B_alg;
    inv_A_theta = inv_A * Theta_mat;
    EinvATheta = E * inv_A_theta;
    if contAngle==1
        N_ini = gradient_simu_N_ini_contact(ind_cont,Pfree_c3,Fc_mc3,angleact,psi_first,contAngle);
        Y = gradient_simu_Y_contact(ind_cont,Rini,Fc_mc3,dPFreec_dp,lp,contAngle,[]);
        EA_1B = E * A_1B + N_ini;
        dOl_dY = EA_1B \ (-EinvATheta + Y);
    else
        Y = gradient_simu_Y_contact(ind_cont,Rini,Fc_mc3,dPFreec_dp,lp,contAngle,dPcq_1_dp);
        EA_1B = (E * A_1B);
        dOl_dY = EA_1B \ (-EinvATheta + Y);
    end
    dF_ddelta_dp = A_1B * dOl_dY + inv_A_theta;
    dFc_dp = dF_ddelta_dp(1:3*length(ind_cont),:);
    ddeltat_dp = dF_ddelta_dp(3*length(ind_cont)+1:end,:);

    dOl_dp = dOl_dY(1:3,:);
    dY_dp = dOl_dY(4:6,:); 
end