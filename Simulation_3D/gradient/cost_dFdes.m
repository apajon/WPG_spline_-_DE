function [der_Kctr_dFdes,der_Kcrot_dFdes] = cost_dFdes(dPc_dFdes,dFc_dFdes,ind_cont,ind_slip,dOl_dFdes,dY_dFdes,ddeltat_dFdes,Cc,P_mc,Wc,Fc_mc,Fc_mc3,PabsOld_mc,fric,contAngle,angleact,displ_first,Rini,R,Pfree_c,dW_Fc_Fdes,Xi,Kc,u_min,u_max,invABwXi,G_inv_A,G_inv_ABw,Fc_mc_hat)
    derG_dFdes = compute_derG_Fdes(dPc_dFdes,dFc_dFdes,ind_cont,ind_slip);
    [derA_dFdes,derWab_i,dR_dFdes]  = compute_derA_Fdes(dY_dFdes,dFc_dFdes,ddeltat_dFdes,Cc,P_mc,Fc_mc,ind_cont,ind_slip,PabsOld_mc,fric,contAngle,angleact,displ_first,Rini,R);
    derBw_dFdes = compute_derBw_Fdes(dFc_dFdes,Fc_mc3,ind_cont,ind_slip,Wc,Pfree_c,derWab_i,dR_dFdes,dW_Fc_Fdes,Fc_mc_hat);
    derXi_dFdes = compute_derXi_Fdes(dOl_dFdes);
    derKc_dFdes = compute_derKc_Fdes(invABwXi,G_inv_A,G_inv_ABw,Xi,derG_dFdes,derA_dFdes,derBw_dFdes,derXi_dFdes);

    %%% For the translational part
    derKc_dFdes_tr = derKc_dFdes(1:3,1:3);
    derKc_dFdes_tr_sym = derKc_dFdes_tr' * Kc(1:3,1:3) + Kc(1:3,1:3)' * derKc_dFdes_tr;
    dlambda1 = u_max'*[derKc_dFdes_tr_sym(1,1) 0 0; 0 0 0; 0 0 0] *u_max;
    dlambda2 = u_max'*[0 derKc_dFdes_tr_sym(1,2) 0; 0 0 0; 0 0 0] *u_max;
    dlambda3 = u_max'*[0 0 derKc_dFdes_tr_sym(1,3); 0 0 0; 0 0 0] *u_max;
    dlambda4 = u_max'*[0 0 0; derKc_dFdes_tr_sym(2,1) 0 0; 0 0 0] *u_max;
    dlambda5 = u_max'*[0 0 0; 0 derKc_dFdes_tr_sym(2,2) 0; 0 0 0] *u_max;
    dlambda6 = u_max'*[0 0 0; 0 0 derKc_dFdes_tr_sym(2,3); 0 0 0] *u_max;
    dlambda7 = u_max'*[0 0 0; 0 0 0; derKc_dFdes_tr_sym(3,1) 0 0] *u_max;
    dlambda8 = u_max'*[0 0 0; 0 0 0; 0 derKc_dFdes_tr_sym(3,2) 0] *u_max;
    dlambda9 = u_max'*[0 0 0; 0 0 0; 0 0 derKc_dFdes_tr_sym(3,3)] *u_max;
    der_Kctr_dFdes = (dlambda1+dlambda2+dlambda3+dlambda4+dlambda5+dlambda6+dlambda7+dlambda8+dlambda9)/(2*sqrt(max(eig(Kc(1:3,1:3)'*Kc(1:3,1:3)))));

    %%% For the rotational part
    derKc_dFdes_rot = derKc_dFdes(4:6,4:6);
    derKc_dFdes_rot_sym = derKc_dFdes_rot' * Kc(4:6,4:6) + Kc(4:6,4:6)' * derKc_dFdes_rot; 
    dlambda1 = u_min'*[derKc_dFdes_rot_sym(1,1) 0 0; 0 0 0; 0 0 0] *u_min;
    dlambda2 = u_min'*[0 derKc_dFdes_rot_sym(1,2) 0; 0 0 0; 0 0 0] *u_min;
    dlambda3 = u_min'*[0 0 derKc_dFdes_rot_sym(1,3); 0 0 0; 0 0 0] *u_min;
    dlambda4 = u_min'*[0 0 0; derKc_dFdes_rot_sym(2,1) 0 0; 0 0 0] *u_min;
    dlambda5 = u_min'*[0 0 0; 0 derKc_dFdes_rot_sym(2,2) 0; 0 0 0] *u_min;
    dlambda6 = u_min'*[0 0 0; 0 0 derKc_dFdes_rot_sym(2,3); 0 0 0] *u_min;
    dlambda7 = u_min'*[0 0 0; 0 0 0; derKc_dFdes_rot_sym(3,1) 0 0] *u_min;
    dlambda8 = u_min'*[0 0 0; 0 0 0; 0 derKc_dFdes_rot_sym(3,2) 0] *u_min;
    dlambda9 = u_min'*[0 0 0; 0 0 0; 0 0 derKc_dFdes_rot_sym(3,3)] *u_min;
    der_Kcrot_dFdes = (dlambda1+dlambda2+dlambda3+dlambda4+dlambda5+dlambda6+dlambda7+dlambda8+dlambda9)/(2*sqrt(min(eig(Kc(4:6,4:6)'*Kc(4:6,4:6)))));

    %%% Gradient of cost function
    %der_cost_dp = der_Kctr_dp - w * der_Kcrot_dp;
end