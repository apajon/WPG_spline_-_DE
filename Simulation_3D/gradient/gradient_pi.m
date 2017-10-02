function [der_Kctr_dpi_tot,der_Kcrot_dpi_tot] = gradient_pi(dPFree_dpi,der_Cs_dpi,sole,fric,displ_first,psi_first,G_tot,A_tot,Bw_tot,Xi_tot,...
                                                            Kcart_tot,u_max_tot,u_min_tot,Ws_tot,PSurf3_tot,Fc_mc3_tot,inv_A_tot,E_tot,...
                                                            A_1B_tot,inv_E_A1B_tot,angleact_tot,ind_cont_tot,ind_slip_tot,Hc_tot)
 der_Kcrot_dpi_tot = 0;
 der_Kctr_dpi_tot = 0;
 for s=1:size(PSurf3_tot,2)
 %for s=1:1
    %%% Init variables %%%
    G = G_tot{s};
    A = A_tot{s};
    Bw = Bw_tot{s};
    Xi = Xi_tot{s};
    Kcart = Kcart_tot(6*s-5:6*s,:);
    u_max = u_max_tot(:,s);
    u_min = u_min_tot(:,s);
    Ws = Ws_tot{s};
    %PSurf3 = PSurf3_tot{s};
    inv_A = inv_A_tot{s};
    E = E_tot{s};
    ind_cont = ind_cont_tot{s};
    Fc_mc3 = Fc_mc3_tot{s};
    A_1B = A_1B_tot{s};
    inv_E_A1B = inv_E_A1B_tot{s};
    angleact = angleact_tot(:,s);
    ind_slip = ind_slip_tot{s};
    Hc = Hc_tot{s};
    theta = angleact(1);
    phi = angleact(2);
    psi = angleact(3);
    R = Rot(theta,phi,psi);
    Rini = Rot(theta,phi,psi_first);
    ind_cont3 = sort(([3*ind_cont-2; 3*ind_cont-1; 3*ind_cont]));
    ind_cont3 = reshape(ind_cont3,3*length(ind_cont),1);      
    der_Cc_dpi = der_Cs_dpi(ind_cont3,ind_cont3);
    dPFrees_dpi = dPFree_dpi(sole.nodesFreeSurf3);
    dPFreec_dpi = dPFrees_dpi(ind_cont3);
    Wc = Ws(ind_cont3,ind_cont3);
    Fc_mc = reshape(Fc_mc3,3,length(ind_cont));
    PSurf = reshape(PSurf3_tot(:,s),3,sole.nFreeSurf);
    P_mc = PSurf(:,ind_cont);
    Pfree = sole.coor(sole.nodesFreeSurf,:)';
    Pfree_s3 = reshape(Pfree,3*sole.nFreeSurf,1);
    Pfree_c = Pfree(:,ind_cont);
    Cs = sole.Cs;
    Cc = Cs(ind_cont3,ind_cont3);
    %%%%%%%%%%%%%%%%%%%%
    if s==1
        PabsOld = Pfree;
        PabsOld_mc = PabsOld(:,ind_cont);
        [Theta_mat_pi,der_HCcH_pi,der_HCcHF_pi] = gradient_simu_Theta_contact_pi(ind_cont,der_Cc_dpi,ind_slip,Fc_mc3,Hc,Rini,dPFreec_dpi,s,[]);           
        Y = gradient_simu_Y_contact_pi(ind_cont,Rini,Fc_mc3,dPFreec_dpi,s,[]);
    else
        PabsOld = reshape(PSurf3_tot(:,s-1),3,sole.nFreeSurf);
        PabsOld_mc = PabsOld(:,ind_cont);
        dPc_dpi_previous = dPs_dpi(ind_cont3);
        [Theta_mat_pi,der_HCcH_pi,der_HCcHF_pi] = gradient_simu_Theta_contact_pi(ind_cont,der_Cc_dpi,ind_slip,Fc_mc3,Hc,Rini,dPFreec_dpi,s,dPc_dpi_previous);           
        Y = gradient_simu_Y_contact_pi(ind_cont,Rini,Fc_mc3,dPFreec_dpi,s,dPc_dpi_previous);          
    end
    inv_A_theta_pi = inv_A * Theta_mat_pi;
    EinvTheta_pi = -E * inv_A_theta_pi;
    inv_E_A1BY_pi = inv_E_A1B * Y;
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
    dW_Fc = Wc * dFc_dpi;
    dPc_dpi =  dPs_dpi(ind_cont3);
    derG_pi = compute_derG_pi(dPc_dpi,dFc_dpi,ind_cont,ind_slip);
    [derA_pi,derWab_i,dR_dpi] = compute_derA_pi(dY_dpi,dFc_dpi,ddeltat_dpi,Cc,P_mc,Fc_mc,ind_cont,ind_slip,PabsOld_mc,fric,s,angleact,displ_first,Rini,der_HCcH_pi,Hc);
    derBw_pi = compute_derBw_pi(dPFreec_dpi,dFc_dpi,Fc_mc3,ind_cont,ind_slip,R,Wc,Pfree_c,derWab_i,dR_dpi,dW_Fc); 
    derXi_pi = compute_derXi_pi(dOl_dpi);
    derKc_pi = compute_derKc_pi(G,A,Bw,Xi,derG_pi,derA_pi,derBw_pi,derXi_pi);

    %%% For the rotational part
    derKc_pi_rot = derKc_pi(4:6,4:6);
    derKc_pi_rot_sym = derKc_pi_rot' * Kcart(4:6,4:6) + Kcart(4:6,4:6)' * derKc_pi_rot; 
    dlambda1 = u_min'*[derKc_pi_rot_sym(1,1) 0 0; 0 0 0; 0 0 0] *u_min;
    dlambda2 = u_min'*[0 derKc_pi_rot_sym(1,2) 0; 0 0 0; 0 0 0] *u_min;
    dlambda3 = u_min'*[0 0 derKc_pi_rot_sym(1,3); 0 0 0; 0 0 0] *u_min;
    dlambda4 = u_min'*[0 0 0; derKc_pi_rot_sym(2,1) 0 0; 0 0 0] *u_min;
    dlambda5 = u_min'*[0 0 0; 0 derKc_pi_rot_sym(2,2) 0; 0 0 0] *u_min;
    dlambda6 = u_min'*[0 0 0; 0 0 derKc_pi_rot_sym(2,3); 0 0 0] *u_min;
    dlambda7 = u_min'*[0 0 0; 0 0 0; derKc_pi_rot_sym(3,1) 0 0] *u_min;
    dlambda8 = u_min'*[0 0 0; 0 0 0; 0 derKc_pi_rot_sym(3,2) 0] *u_min;
    dlambda9 = u_min'*[0 0 0; 0 0 0; 0 0 derKc_pi_rot_sym(3,3)] *u_min;
    der_Kcrot_dpi = (dlambda1+dlambda2+dlambda3+dlambda4+dlambda5+dlambda6+dlambda7+dlambda8+dlambda9)/(2*sqrt(min(eig(Kcart(4:6,4:6)'*Kcart(4:6,4:6)))));
    if s<(size(PSurf3_tot,2)/10)
        %%% For the translational part
        derKc_pi_tr = derKc_pi(1:3,1:3);
        derKc_pi_tr_sym = derKc_pi_tr' * Kcart(1:3,1:3) + Kcart(1:3,1:3)' * derKc_pi_tr;
        dlambda1 = u_max'*[derKc_pi_tr_sym(1,1) 0 0; 0 0 0; 0 0 0] *u_max;
        dlambda2 = u_max'*[0 derKc_pi_tr_sym(1,2) 0; 0 0 0; 0 0 0] *u_max;
        dlambda3 = u_max'*[0 0 derKc_pi_tr_sym(1,3); 0 0 0; 0 0 0] *u_max;
        dlambda4 = u_max'*[0 0 0; derKc_pi_tr_sym(2,1) 0 0; 0 0 0] *u_max;
        dlambda5 = u_max'*[0 0 0; 0 derKc_pi_tr_sym(2,2) 0; 0 0 0] *u_max;
        dlambda6 = u_max'*[0 0 0; 0 0 derKc_pi_tr_sym(2,3); 0 0 0] *u_max;
        dlambda7 = u_max'*[0 0 0; 0 0 0; derKc_pi_tr_sym(3,1) 0 0] *u_max;
        dlambda8 = u_max'*[0 0 0; 0 0 0; 0 derKc_pi_tr_sym(3,2) 0] *u_max;
        dlambda9 = u_max'*[0 0 0; 0 0 0; 0 0 derKc_pi_tr_sym(3,3)] *u_max;
        der_Kctr_dpi = (dlambda1+dlambda2+dlambda3+dlambda4+dlambda5+dlambda6+dlambda7+dlambda8+dlambda9)/(2*sqrt(max(eig(Kcart(1:3,1:3)'*Kcart(1:3,1:3)))));
        der_Kctr_dpi_tot = der_Kctr_dpi_tot + der_Kctr_dpi;
    end
    der_Kcrot_dpi_tot = der_Kcrot_dpi_tot + der_Kcrot_dpi;
end
end