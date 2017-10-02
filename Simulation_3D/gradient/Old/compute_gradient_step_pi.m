function [der_Kctr_dpi,der_Kcrot_dpi,dPs_dpi] = compute_gradient_step_pi(sole,G,A,Bw,Xi,Kcart,u_max,u_min,ind_cont3,Ws,Wc,Cc,P_mc,Fc_mc,PabsOld_mc,der_Cs_dpi,A_1B,EinvTheta_pi,inv_E_A1B,inv_E_A1BY_pi,inv_A_theta_pi,Pfree_s3,Pfree_c,angleact,fric,Fc_mc3,Cs,ind_cont,ind_slip,displ_first,contAngle,dPFrees_dpi,dPFreec_dpi,R,Rini,der_HCcH_pi,Hc,w)
    dOl_dY = inv_E_A1B * EinvTheta_pi + inv_E_A1BY_pi;
    dF_ddelta_dp = A_1B * dOl_dY + inv_A_theta_pi;
    dFc_dpi = dF_ddelta_dp(1:3*length(ind_cont));
    ddeltat_dpi = dF_ddelta_dp(3*length(ind_cont)+1:end);

    dOl_dpi = dOl_dY(1:3);
    dY_dpi = dOl_dY(4:6);     
   
    dPs_dpi = zeros(3*sole.nFreeSurf,1);
    [dR_dtheta,dR_dphi,dR_dpsi] = dRot(angleact);
    for j=1:sole.nFreeSurf
        P = Pfree_s3((3*j-2):3*j);
        btheta = dR_dtheta*P;
        bphi = dR_dphi*P;
        bpsi = dR_dpsi*P;
        W_Fc = zeros(3,1);
        RdCRt = zeros(3,1);
        for k=1:length(ind_cont)
            F = Fc_mc3(3*k-2:3*k);
            C = Cs(3*j-2:3*j,3*ind_cont(k)-2:3*ind_cont(k));
            RtF = R'*F;
            CRtF = C*RtF;
            btheta = btheta + dR_dtheta * CRtF + R * (C * (dR_dtheta' * F));
            bphi = bphi + dR_dphi * CRtF + R * (C * (dR_dphi' * F));
            bpsi = bpsi + dR_dpsi  * CRtF + R * (C * (dR_dpsi' * F));
            W_Fc = W_Fc + Ws(3*j-2:3*j,3*ind_cont(k)-2:3*ind_cont(k)) * dFc_dpi((3*k-2):3*k);
            RdCRt = RdCRt + R * (der_Cs_dpi(3*j-2:3*j,3*ind_cont(k)-2:3*ind_cont(k)) * RtF);                  
        end
        dPs_dpi(3*j-2:3*j) = dOl_dpi + [btheta bphi bpsi] * dY_dpi + R * dPFrees_dpi(3*j-2:3*j) + W_Fc + RdCRt;
    end  
    %[dPs_dp,dW_Fc] = derPs_dp(sole,angleact,ind_cont,Wc,Pfree_s3,dOl_dp,dY_dp,dFc_dp,lp,dPFrees_dp,B12out,der_HCcHF);
    %[dPs_dp] = derPs_dp(sole,angleact,ind_cont,Ws,Pfree_s3,dOl_dp,dY_dp,dFc_dp,lp,dPFrees_dp,Cs,sole.der_Cs_dp,Fc_mc3);
    %C = setdiff(A,B);
    dW_Fc = Wc * dFc_dpi;
    dPc_dpi =  dPs_dpi(ind_cont3);
    %der_cost_dp = zeros(lp,1); 
%     der_Kctr_dp = zeros(lp,1);
%     der_Kcrot_dp = zeros(lp,1);
%     for i=1:lp
%         [der_Kctr_dp(i,1),der_Kcrot_dp(i,1)] = cost_dp(dPc_dp(:,i),dFc_dp(:,i),ind_cont,ind_slip,dOl_dp(:,i),dY_dp(:,i),ddeltat_dp(:,i),Cc,P_mc,Wc,Fc_mc,Fc_mc3,PabsOld_mc,fric,contAngle,angleact,displ_first,Rini,der_HCcH{i},Hc,dPFreec_dp(:,i),R,Pfree_c,dW_Fc(:,i),G,A,Bw,Xi,Kcart,w,u_min,u_max);
%     end
    
    derG_p_i = derG_pi(dPc_dpi,dFc_dpi,ind_cont,ind_slip);
    [derA_p_i,derWab_i,dR_dp_i] = derA_pi(dY_dpi,dFc_dpi,ddeltat_dpi,Cc,P_mc,Fc_mc,ind_cont,ind_slip,PabsOld_mc,fric,contAngle,angleact,displ_first,Rini,der_HCcH_pi,Hc);
    derBw_p_i = derBw_pi(dPFreec_dpi,dFc_dpi,Fc_mc3,ind_cont,ind_slip,R,Wc,Pfree_c,derWab_i,dR_dp_i,dW_Fc); 
    derXi_p_i = derXi_pi(dOl_dpi);
    derKc_p_i = derKc_pi(G,A,Bw,Xi,derG_p_i,derA_p_i,derBw_p_i,derXi_p_i);

    %%% For the translational part
    derKc_p_i_tr = derKc_p_i(1:3,1:3);
    derKc_p_i_tr_sym = derKc_p_i_tr' * Kcart(1:3,1:3) + Kcart(1:3,1:3)' * derKc_p_i_tr;
    dlambda1 = u_max'*[derKc_p_i_tr_sym(1,1) 0 0; 0 0 0; 0 0 0] *u_max;
    dlambda2 = u_max'*[0 derKc_p_i_tr_sym(1,2) 0; 0 0 0; 0 0 0] *u_max;
    dlambda3 = u_max'*[0 0 derKc_p_i_tr_sym(1,3); 0 0 0; 0 0 0] *u_max;
    dlambda4 = u_max'*[0 0 0; derKc_p_i_tr_sym(2,1) 0 0; 0 0 0] *u_max;
    dlambda5 = u_max'*[0 0 0; 0 derKc_p_i_tr_sym(2,2) 0; 0 0 0] *u_max;
    dlambda6 = u_max'*[0 0 0; 0 0 derKc_p_i_tr_sym(2,3); 0 0 0] *u_max;
    dlambda7 = u_max'*[0 0 0; 0 0 0; derKc_p_i_tr_sym(3,1) 0 0] *u_max;
    dlambda8 = u_max'*[0 0 0; 0 0 0; 0 derKc_p_i_tr_sym(3,2) 0] *u_max;
    dlambda9 = u_max'*[0 0 0; 0 0 0; 0 0 derKc_p_i_tr_sym(3,3)] *u_max;
    der_Kctr_dpi = (dlambda1+dlambda2+dlambda3+dlambda4+dlambda5+dlambda6+dlambda7+dlambda8+dlambda9)/(2*sqrt(max(eig(Kcart(1:3,1:3)'*Kcart(1:3,1:3)))));

    %%% For the rotational part
    derKc_p_i_rot = derKc_p_i(4:6,4:6);
    derKc_p_i_rot_sym = derKc_p_i_rot' * Kcart(4:6,4:6) + Kcart(4:6,4:6)' * derKc_p_i_rot; 
    dlambda1 = u_min'*[derKc_p_i_rot_sym(1,1) 0 0; 0 0 0; 0 0 0] *u_min;
    dlambda2 = u_min'*[0 derKc_p_i_rot_sym(1,2) 0; 0 0 0; 0 0 0] *u_min;
    dlambda3 = u_min'*[0 0 derKc_p_i_rot_sym(1,3); 0 0 0; 0 0 0] *u_min;
    dlambda4 = u_min'*[0 0 0; derKc_p_i_rot_sym(2,1) 0 0; 0 0 0] *u_min;
    dlambda5 = u_min'*[0 0 0; 0 derKc_p_i_rot_sym(2,2) 0; 0 0 0] *u_min;
    dlambda6 = u_min'*[0 0 0; 0 0 derKc_p_i_rot_sym(2,3); 0 0 0] *u_min;
    dlambda7 = u_min'*[0 0 0; 0 0 0; derKc_p_i_rot_sym(3,1) 0 0] *u_min;
    dlambda8 = u_min'*[0 0 0; 0 0 0; 0 derKc_p_i_rot_sym(3,2) 0] *u_min;
    dlambda9 = u_min'*[0 0 0; 0 0 0; 0 0 derKc_p_i_rot_sym(3,3)] *u_min;
    der_Kcrot_dpi = (dlambda1+dlambda2+dlambda3+dlambda4+dlambda5+dlambda6+dlambda7+dlambda8+dlambda9)/(2*sqrt(min(eig(Kcart(4:6,4:6)'*Kcart(4:6,4:6)))));

    
end