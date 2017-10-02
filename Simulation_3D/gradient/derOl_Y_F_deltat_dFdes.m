function  [dOl_dFdes,dY_dFdes,dFc_dFdes,ddeltat_dFdes] = derOl_Y_F_deltat_dFdes(ind_cont,ind_slip,E,inv_A,EA_1B,A_1B,contAngle,dPc_q_1_dFdes,Fc_mc3,ind)
    Y = gradient_simu_Y_contact_dFdes(dPc_q_1_dFdes,contAngle,ind_cont,Fc_mc3,ind);
    if contAngle==ind
        dOl_dY = EA_1B \ Y;
        dF_ddelta_dFdes = A_1B * dOl_dY;
    else
        Theta_mat = gradient_simu_Theta_contact_dFdes(ind_slip,dPc_q_1_dFdes);
        dOl_dY = EA_1B \ (- E * inv_A * Theta_mat + Y);
        dF_ddelta_dFdes = A_1B * dOl_dY + inv_A * Theta_mat;
    end
    dFc_dFdes = dF_ddelta_dFdes(1:3*length(ind_cont),:);
    ddeltat_dFdes = dF_ddelta_dFdes(3*length(ind_cont)+1:end,:);

    dOl_dFdes = dOl_dY(1:3,:);
    dY_dFdes = dOl_dY(4:6,:);
end