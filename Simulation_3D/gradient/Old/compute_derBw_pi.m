function derBw_p_i = compute_derBw_pi(dFc_dp_i,Fc_mc3,ind_cont,ind_slip,Wc,Pfree_c,derWab_i,dR_dp_i,dW_Fci,Fc_mc_hat,RdPFreec_dp)
    %dFc_mc_hat_dp_i = zeros(3*length(ind_cont),3);
%     for k=1:length(ind_cont)        
%         dFc_mc_hat_dp_i(3*k-2:3*k,:) = skew(dFc_dp_i(3*k-2:3*k));                         
%     end
    dFc_mc_hat_dp_i = multiSkew(dFc_dp_i);
    %derBw_p_i12 = zeros(3*length(ind_cont),3);
    dWFh = derWab_i * Fc_mc_hat;
    WdF = Wc * dFc_mc_hat_dp_i;
    dWf3 = derWab_i * Fc_mc3;
%     for k=1:length(ind_cont)
%         b = dWFh(3*k-2:3*k,:);
%         c = WdF(3*k-2:3*k,:);
%         d = dWf3(3*k-2:3*k) + dW_Fci(3*k-2:3*k);
%         a = dR_dp_i * Pfree_c(:,k) + RdPFreec_dp(3*k-2:3*k);
%         ad_skew = skew(a+d);
%         derBw_p_i12(3*k-2:3*k,:) = -ad_skew + b + c;       
%     end
    a = reshape(dR_dp_i*Pfree_c,3*length(ind_cont),1) + RdPFreec_dp;
    ad = a + dWf3 + dW_Fci;
    derBw_p_i12 = dWFh + WdF - multiSkew(ad);
    derBw_p_i = [zeros(3*length(ind_cont),3) derBw_p_i12;zeros(2*length(ind_slip),3) zeros(2*length(ind_slip),3)];    
end