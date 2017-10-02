function  [dOl_dZdes,dY_dZdes,dFc_dZdes,ddeltat_dZdes] = derOl_Y_F_deltat_dZdes(ind_cont,ind_slip,E,inv_A,EA_1B,A_1B,contAngle,dPc_q_1_dZdes,Fc_mc3,ind)
    Y = gradient_simu_Y_contact_dZdes(dPc_q_1_dZdes,contAngle,ind_cont,Fc_mc3,ind);
    if contAngle==ind
        dOl_dY = EA_1B \ Y;
        dF_ddelta_dZdes = A_1B * dOl_dY;
    else
        Theta_mat = gradient_simu_Theta_contact_dZdes(ind_slip,dPc_q_1_dZdes);
        dOl_dY = EA_1B \ (- E * inv_A * Theta_mat + Y);
        dF_ddelta_dZdes = A_1B * dOl_dY + inv_A * Theta_mat;
    end
    dFc_dZdes = dF_ddelta_dZdes(1:3*length(ind_cont),:);
    ddeltat_dZdes = dF_ddelta_dZdes(3*length(ind_cont)+1:end,:);

    dOl_dZdes = dOl_dY(1:3,:);
    dY_dZdes = dOl_dY(4:6,:);
end