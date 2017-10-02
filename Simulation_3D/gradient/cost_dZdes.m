function [der_Kctr_dZdes,der_Kcrot_dZdes] = cost_dZdes(index_Zdes,dPc_dZdes,dFc_dZdes,ind_cont,ind_slip,dOl_dZdes,dY_dZdes,ddeltat_dZdes,Cc,P_mc,Wc,Fc_mc,Fc_mc3,PabsOld_mc,fric,contAngle,angleact,displ_first,Rini,R,Pfree_c,dW_Fc_Zdes,Xi,Kc,u_min,u_max,invABwXi,G_inv_A,G_inv_ABw,Fc_mc_hat,ind)
    derG_dZdes = compute_derG_Zdes(dPc_dZdes,dFc_dZdes,ind_cont,ind_slip,index_Zdes,ind,contAngle);
    [derA_dZdes,derWab_i,dR_dZdes]  = compute_derA_Zdes(dY_dZdes,dFc_dZdes,ddeltat_dZdes,Cc,P_mc,Fc_mc,ind_cont,ind_slip,PabsOld_mc,fric,contAngle,angleact,displ_first,Rini,R);
    derBw_dZdes = compute_derBw_Zdes(dFc_dZdes,Fc_mc3,ind_cont,ind_slip,Wc,Pfree_c,derWab_i,dR_dZdes,dW_Fc_Zdes,Fc_mc_hat);
    derXi_dZdes = compute_derXi_Zdes(dOl_dZdes,index_Zdes,ind,contAngle);
    derKc_dZdes = compute_derKc_Zdes(invABwXi,G_inv_A,G_inv_ABw,Xi,derG_dZdes,derA_dZdes,derBw_dZdes,derXi_dZdes);

    %%% For the translational part
    derKc_dZdes_tr = derKc_dZdes(1:3,1:3);
    derKc_dZdes_tr_sym = derKc_dZdes_tr' * Kc(1:3,1:3) + Kc(1:3,1:3)' * derKc_dZdes_tr;
    dlambda1 = u_max'*[derKc_dZdes_tr_sym(1,1) 0 0; 0 0 0; 0 0 0] *u_max;
    dlambda2 = u_max'*[0 derKc_dZdes_tr_sym(1,2) 0; 0 0 0; 0 0 0] *u_max;
    dlambda3 = u_max'*[0 0 derKc_dZdes_tr_sym(1,3); 0 0 0; 0 0 0] *u_max;
    dlambda4 = u_max'*[0 0 0; derKc_dZdes_tr_sym(2,1) 0 0; 0 0 0] *u_max;
    dlambda5 = u_max'*[0 0 0; 0 derKc_dZdes_tr_sym(2,2) 0; 0 0 0] *u_max;
    dlambda6 = u_max'*[0 0 0; 0 0 derKc_dZdes_tr_sym(2,3); 0 0 0] *u_max;
    dlambda7 = u_max'*[0 0 0; 0 0 0; derKc_dZdes_tr_sym(3,1) 0 0] *u_max;
    dlambda8 = u_max'*[0 0 0; 0 0 0; 0 derKc_dZdes_tr_sym(3,2) 0] *u_max;
    dlambda9 = u_max'*[0 0 0; 0 0 0; 0 0 derKc_dZdes_tr_sym(3,3)] *u_max;
    der_Kctr_dZdes = (dlambda1+dlambda2+dlambda3+dlambda4+dlambda5+dlambda6+dlambda7+dlambda8+dlambda9)/(2*sqrt(max(eig(Kc(1:3,1:3)'*Kc(1:3,1:3)))));

    %%% For the rotational part
    derKc_dZdes_rot = derKc_dZdes(4:6,4:6);
    derKc_dZdes_rot_sym = derKc_dZdes_rot' * Kc(4:6,4:6) + Kc(4:6,4:6)' * derKc_dZdes_rot; 
    dlambda1 = u_min'*[derKc_dZdes_rot_sym(1,1) 0 0; 0 0 0; 0 0 0] *u_min;
    dlambda2 = u_min'*[0 derKc_dZdes_rot_sym(1,2) 0; 0 0 0; 0 0 0] *u_min;
    dlambda3 = u_min'*[0 0 derKc_dZdes_rot_sym(1,3); 0 0 0; 0 0 0] *u_min;
    dlambda4 = u_min'*[0 0 0; derKc_dZdes_rot_sym(2,1) 0 0; 0 0 0] *u_min;
    dlambda5 = u_min'*[0 0 0; 0 derKc_dZdes_rot_sym(2,2) 0; 0 0 0] *u_min;
    dlambda6 = u_min'*[0 0 0; 0 0 derKc_dZdes_rot_sym(2,3); 0 0 0] *u_min;
    dlambda7 = u_min'*[0 0 0; 0 0 0; derKc_dZdes_rot_sym(3,1) 0 0] *u_min;
    dlambda8 = u_min'*[0 0 0; 0 0 0; 0 derKc_dZdes_rot_sym(3,2) 0] *u_min;
    dlambda9 = u_min'*[0 0 0; 0 0 0; 0 0 derKc_dZdes_rot_sym(3,3)] *u_min;
    der_Kcrot_dZdes = (dlambda1+dlambda2+dlambda3+dlambda4+dlambda5+dlambda6+dlambda7+dlambda8+dlambda9)/(2*sqrt(min(eig(Kc(4:6,4:6)'*Kc(4:6,4:6)))));

    %%% Gradient of cost function
    %der_cost_dp = der_Kctr_dp - w * der_Kcrot_dp;
end