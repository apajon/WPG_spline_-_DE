function Theta_mat = gradient_simu_Theta_contact_dFdes(ind_slip,dPc_q_1_dFdes)
    dPc_q_1_dFdes(3:3:end,:) = 0;

    Theta_mat = [-dPc_q_1_dFdes;zeros(2*length(ind_slip),3)];
end