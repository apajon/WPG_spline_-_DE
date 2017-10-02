function derBw_Fdes = compute_derBw_Fdes(dFc_dFdes,Fc_mc3,ind_cont,ind_slip,Wc,Pfree_c,derWab_i,dR_dFdes,dW_Fc_Fdes,Fc_mc_hat)
    dFc_mc_hat_dFdes = multiSkew(dFc_dFdes);
    dWFh = derWab_i * Fc_mc_hat;
    WdF = Wc * dFc_mc_hat_dFdes;
    dWf3 = derWab_i * Fc_mc3;
    a = reshape(dR_dFdes*Pfree_c,3*length(ind_cont),1);
    ad = a + dWf3 + dW_Fc_Fdes;
    derBw_p_i12 = dWFh + WdF - multiSkew(ad);
    derBw_Fdes = [zeros(3*length(ind_cont),3) derBw_p_i12;zeros(2*length(ind_slip),3) zeros(2*length(ind_slip),3)];    
end