function  [dFc_dFdes,ddeltat_dFdes] = derOl_Y_F_deltat_dFdes1(ind_cont,ind_slip,E,inv_A,EA_1B,A_1B,contAngle,dPc_q_1_dFdes,Fc_mc3)
    dPc_q_1_dFdes(3:3:end,:) = 0;
    Theta_mat = gradient_simu_Theta_contact_dFdes(ind_slip,dPc_q_1_dFdes);
    dF_ddelta_dFdes = inv_A * Theta_mat;
    dFc_dFdes = dF_ddelta_dFdes(1:3*length(ind_cont),:);
    ddeltat_dFdes = dF_ddelta_dFdes(3*length(ind_cont)+1:end,:);
end