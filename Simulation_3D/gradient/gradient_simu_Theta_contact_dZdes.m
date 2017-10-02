function Theta_mat = gradient_simu_Theta_contact_dZdes(ind_slip,dPc_q_1_dZdes)
    dPc_q_1_dZdes(3:3:end,:) = 0;

    Theta_mat = [-dPc_q_1_dZdes;zeros(2*length(ind_slip),2)];
end