function [der_Kctr_dp,der_Kcrot_dp,derKc_p_i] = cost_dp(dPc_dp,dFc_dp,ind_cont,ind_slip,dOl_dp,dY_dp,ddeltat_dp,Cc,P_mc,Wc,Fc_mc,Fc_mc3,PabsOld_mc,fric,contAngle,angleact,displ_first,Rini,der_HCcH,R,Pfree_c,dW_Fc,Xi,Kc,u_min,u_max,invABwXi,G_inv_A,G_inv_ABw,Fc_mc_hat,RdPFreec_dp)
    derG_p_i = compute_derG_pi(dPc_dp,dFc_dp,ind_cont,ind_slip);
    [derA_p_i,derWab_i,dR_dp_i] = compute_derA_pi(dY_dp,dFc_dp,ddeltat_dp,Cc,P_mc,Fc_mc,ind_cont,ind_slip,PabsOld_mc,fric,contAngle,angleact,displ_first,Rini,der_HCcH,R);
    derBw_p_i = compute_derBw_pi(dFc_dp,Fc_mc3,ind_cont,ind_slip,Wc,Pfree_c,derWab_i,dR_dp_i,dW_Fc,Fc_mc_hat,RdPFreec_dp); 
    derXi_p_i = compute_derXi_pi(dOl_dp);
    derKc_p_i = compute_derKc_pi(invABwXi,G_inv_A,G_inv_ABw,Xi,derG_p_i,derA_p_i,derBw_p_i,derXi_p_i);

    %%% For the translational part
    derKc_p_i_tr = derKc_p_i(1:3,1:3);
    derKc_p_i_tr_sym = derKc_p_i_tr' * Kc(1:3,1:3) + Kc(1:3,1:3)' * derKc_p_i_tr;
    dlambda1 = u_max'*[derKc_p_i_tr_sym(1,1) 0 0; 0 0 0; 0 0 0] *u_max;
    dlambda2 = u_max'*[0 derKc_p_i_tr_sym(1,2) 0; 0 0 0; 0 0 0] *u_max;
    dlambda3 = u_max'*[0 0 derKc_p_i_tr_sym(1,3); 0 0 0; 0 0 0] *u_max;
    dlambda4 = u_max'*[0 0 0; derKc_p_i_tr_sym(2,1) 0 0; 0 0 0] *u_max;
    dlambda5 = u_max'*[0 0 0; 0 derKc_p_i_tr_sym(2,2) 0; 0 0 0] *u_max;
    dlambda6 = u_max'*[0 0 0; 0 0 derKc_p_i_tr_sym(2,3); 0 0 0] *u_max;
    dlambda7 = u_max'*[0 0 0; 0 0 0; derKc_p_i_tr_sym(3,1) 0 0] *u_max;
    dlambda8 = u_max'*[0 0 0; 0 0 0; 0 derKc_p_i_tr_sym(3,2) 0] *u_max;
    dlambda9 = u_max'*[0 0 0; 0 0 0; 0 0 derKc_p_i_tr_sym(3,3)] *u_max;
    der_Kctr_dp = (dlambda1+dlambda2+dlambda3+dlambda4+dlambda5+dlambda6+dlambda7+dlambda8+dlambda9)/(2*sqrt(max(eig(Kc(1:3,1:3)'*Kc(1:3,1:3)))));

    %%% For the rotational part
    derKc_p_i_rot = derKc_p_i(4:6,4:6);
    derKc_p_i_rot_sym = derKc_p_i_rot' * Kc(4:6,4:6) + Kc(4:6,4:6)' * derKc_p_i_rot; 
    dlambda1 = u_min'*[derKc_p_i_rot_sym(1,1) 0 0; 0 0 0; 0 0 0] *u_min;
    dlambda2 = u_min'*[0 derKc_p_i_rot_sym(1,2) 0; 0 0 0; 0 0 0] *u_min;
    dlambda3 = u_min'*[0 0 derKc_p_i_rot_sym(1,3); 0 0 0; 0 0 0] *u_min;
    dlambda4 = u_min'*[0 0 0; derKc_p_i_rot_sym(2,1) 0 0; 0 0 0] *u_min;
    dlambda5 = u_min'*[0 0 0; 0 derKc_p_i_rot_sym(2,2) 0; 0 0 0] *u_min;
    dlambda6 = u_min'*[0 0 0; 0 0 derKc_p_i_rot_sym(2,3); 0 0 0] *u_min;
    dlambda7 = u_min'*[0 0 0; 0 0 0; derKc_p_i_rot_sym(3,1) 0 0] *u_min;
    dlambda8 = u_min'*[0 0 0; 0 0 0; 0 derKc_p_i_rot_sym(3,2) 0] *u_min;
    dlambda9 = u_min'*[0 0 0; 0 0 0; 0 0 derKc_p_i_rot_sym(3,3)] *u_min;
    der_Kcrot_dp = (dlambda1+dlambda2+dlambda3+dlambda4+dlambda5+dlambda6+dlambda7+dlambda8+dlambda9)/(2*sqrt(min(eig(Kc(4:6,4:6)'*Kc(4:6,4:6)))));

    %%% Gradient of cost function
    %der_cost_dp = der_Kctr_dp - w * der_Kcrot_dp;
end