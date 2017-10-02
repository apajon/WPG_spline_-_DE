function derBw_Zdes = compute_derBw_Zdes(dFc_dZdes,Fc_mc3,ind_cont,ind_slip,Wc,Pfree_c,derWab_i,dR_dZdes,dW_Fc_Zdes,Fc_mc_hat)
    dFc_mc_hat_dZdes = multiSkew(dFc_dZdes);
    dWFh = derWab_i * Fc_mc_hat;
    WdF = Wc * dFc_mc_hat_dZdes;
    dWf3 = derWab_i * Fc_mc3;
    a = reshape(dR_dZdes*Pfree_c,3*length(ind_cont),1);
    ad = a + dWf3 + dW_Fc_Zdes;
    derBw_p_i12 = dWFh + WdF - multiSkew(ad);
    derBw_Zdes = [zeros(3*length(ind_cont),3) derBw_p_i12;zeros(2*length(ind_slip),3) zeros(2*length(ind_slip),3)];    
end