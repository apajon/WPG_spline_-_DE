function [der_Kctr_dp,der_Kcrot_dp,dPs_dp,der_Kctr_dFdes,der_Kcrot_dFdes,dPs_dFdes,der_Kctr_dZdes,der_Kcrot_dZdes,dPs_dZdes] = compute_gradient_step(sole,lp,Kcart,displ,angleact,Zdes,fric,Pfree_s,PSurf3,PabsOld,Fc_mc3,Cs,dPFree_dp,ind_cont,ind_slip,displ_first,psi_first,contAngle,w,B_alg,dPs_dp,dPs_dFdes,dPs_dZdes)
    D = 3;
    ind_cont3 = sort(([D*ind_cont-2; D*ind_cont-1; D*ind_cont]));
    ind_cont3 = reshape(ind_cont3,D*length(ind_cont),1);    
    theta = angleact(1);
    phi = angleact(2);
    psi = angleact(3);
    R = Rot(theta,phi,psi);
    Rini = Rot(theta,phi,psi_first);
    Pfree_c = Pfree_s(:,ind_cont);
    Pfree_c3 = reshape(Pfree_c,3*length(ind_cont),1);
    Pfree_s3 = reshape(Pfree_s,3*sole.nFreeSurf,1);
    Fc_mc = reshape(Fc_mc3,3,length(ind_cont));
    P_mc3 = PSurf3(ind_cont3);
    P_mc = reshape(P_mc3,3,length(ind_cont)); 
%     Ws = zeros(3*sole.nFreeSurf,3*sole.nFreeSurf);
%     for j = 1:sole.nFreeSurf
%        for k = 1:sole.nFreeSurf
%            Ws(3*j-2:3*j,(3*k)-2:3*k) = R * Cs(3*j-2:3*j,3*k-2:3*k) * R';
%        end
%     end
    nf = sole.nFreeSurf;
    tmp = reshape(R*reshape(Cs,3,3*nf*nf),3*nf,3*nf);
    tmp2 = tmp';
    Ws = reshape(R*reshape(tmp2,3,3*nf*nf),3*nf,3*nf);
    %Ws = reshape(R*reshape(reshape(Cs,3*nf*nf,3)*R',3,3*nf*nf),3*nf,3*nf);
    Cc = Cs(ind_cont3,ind_cont3);
    Wc = Ws(ind_cont3,ind_cont3);
    PabsOld_mc = PabsOld(:,ind_cont);
    der_Cc_dp = cell(lp,1);
    for i = 1:lp
        der_Cc_dp{i} = sole.der_Cs_dp{i}(ind_cont3,ind_cont3);
    end
    dPFrees_dp = dPFree_dp(sole.nodesFreeSurf3,:);
    dPFreec_dp = dPFrees_dp(ind_cont3,:);
    %%% Compute Kcart     
    A = gradient_simu_A_contact(P_mc,Fc_mc,ind_cont,Wc,ind_slip,PabsOld_mc,fric,contAngle,Rini,displ_first);
    G = computeG(ind_cont,ind_slip,Fc_mc,P_mc,Zdes);
    [Bw,Fc_mc_hat] = computeBw(Fc_mc3,ind_cont,ind_slip,R,Wc,Pfree_c);
    Xi = computeXi(displ,Zdes);
    %Kcart = computeKc(G,A,Bw,Xi);    
    %%% Compute gradient of Kcart
    if contAngle==1
        dPcq_1_dp = zeros(3*length(ind_cont),lp);
    else
        dPcq_1_dp = dPs_dp(ind_cont3,:);
    end
    [dOl_dp,dY_dp,dFc_dp,ddeltat_dp,der_HCcH,inv_A,A_1B,EA_1B,E] = derOl_Y_F_deltat_dp(A,Zdes,angleact,Cc,der_Cc_dp,ind_cont,ind_slip,Fc_mc,Fc_mc3,P_mc,Pfree_c3,contAngle,psi_first,lp,dPFreec_dp,B_alg,dPcq_1_dp);
        
    [dPs_dp,RdPFrees_dp,b] = derPs_dp(sole,angleact,ind_cont,ind_cont3,Ws,Pfree_s3,dOl_dp,dY_dp,dFc_dp,lp,dPFrees_dp,Cs,sole.der_Cs_dp,Fc_mc3);
      
    RdPFreec_dp = RdPFrees_dp(ind_cont3);

    dW_Fc = Wc * dFc_dp;    
    dPc_dp =  dPs_dp(ind_cont3,:);

    Kc_squared_tr = Kcart(1:3,1:3)'*Kcart(1:3,1:3);
    [u,DE] = eig(Kc_squared_tr);
    [lambda_max,I_lamba_max] = max(diag(DE));
    u_max = u(:,I_lamba_max);
    Kc_squared_rot = Kcart(4:6,4:6)'*Kcart(4:6,4:6);
    [u,DE] = eig(Kc_squared_rot);
    [lambda_min,I_lamba_min] = min(diag(DE));
    u_min = u(:,I_lamba_min); 
    invABwXi = inv_A * Bw * Xi;
    G_inv_A = G * inv_A;
    G_inv_ABw = G * inv_A * Bw;
    der_Kctr_dp = zeros(lp,1);
    der_Kcrot_dp = zeros(lp,1);
    derKc_p = cell(lp,1);
    for i=1:lp
        [der_Kctr_dp(i),der_Kcrot_dp(i),derKc_p{i}] = cost_dp(dPc_dp(:,i),dFc_dp(:,i),ind_cont,ind_slip,dOl_dp(:,i),dY_dp(:,i),ddeltat_dp(:,i),Cc,P_mc,Wc,Fc_mc,Fc_mc3,PabsOld_mc,fric,contAngle,angleact,displ_first,Rini,der_HCcH{i},R,Pfree_c,dW_Fc(:,i),Xi,Kcart,u_min,u_max,invABwXi,G_inv_A,G_inv_ABw,Fc_mc_hat,RdPFreec_dp);
    end
   
    %%% compute Fdes Gradient
    der_Kctr_dFdes = cell(length(1:contAngle),1);
    der_Kcrot_dFdes = cell(length(1:contAngle),1);
    for j=1:contAngle
        if contAngle==1
            dPc_q_1_dFdes = [];            
        else
            if j < contAngle
                dPc_q_1_dFdes = dPs_dFdes{j}(ind_cont3,:);
            else
                dPc_q_1_dFdes = [];
            end
        end
        [dOl_dFdes,dY_dFdes,dFc_dFdes,ddeltat_dFdes] = derOl_Y_F_deltat_dFdes(ind_cont,ind_slip,E,inv_A,EA_1B,A_1B,contAngle,dPc_q_1_dFdes,Fc_mc3,j);
        dPs_dFdes{j} = derPs_dFdes(sole,ind_cont3,Ws,dOl_dFdes,dY_dFdes,dFc_dFdes,b);
        dW_Fc_Fdes = Wc * dFc_dFdes;
        dPc_dFdes =  dPs_dFdes{j}(ind_cont3,:);
        for i=1:3
            [der_Kctr_dFdes{j}(i,1),der_Kcrot_dFdes{j}(i,1)] = cost_dFdes(dPc_dFdes(:,i),dFc_dFdes(:,i),ind_cont,ind_slip,dOl_dFdes(:,i),dY_dFdes(:,i),ddeltat_dFdes(:,i),Cc,P_mc,Wc,Fc_mc,Fc_mc3,PabsOld_mc,fric,contAngle,angleact,displ_first,Rini,R,Pfree_c,dW_Fc_Fdes(:,i),Xi,Kcart,u_min,u_max,invABwXi,G_inv_A,G_inv_ABw,Fc_mc_hat);   
        end        
    end
    
    %%% compute Zdes Gradient
    der_Kctr_dZdes = cell(length(1:contAngle),1);
    der_Kcrot_dZdes = cell(length(1:contAngle),1);    
    for j=1:contAngle
        if contAngle==1
            dPc_q_1_dZdes = [];            
        else
            if j < contAngle
                dPc_q_1_dZdes = dPs_dZdes{j}(ind_cont3,:);
            else
                dPc_q_1_dZdes = [];
            end
        end
        %[dOl_dFdes,dY_dFdes,dFc_dFdes,ddeltat_dFdes] = derOl_Y_F_deltat_dFdes(ind_cont,ind_slip,E,inv_A,EA_1B,A_1B,contAngle,dPc_q_1_dFdes,Fc_mc3,j);
        [dOl_dZdes,dY_dZdes,dFc_dZdes,ddeltat_dZdes] = derOl_Y_F_deltat_dZdes(ind_cont,ind_slip,E,inv_A,EA_1B,A_1B,contAngle,dPc_q_1_dZdes,Fc_mc3,j);
        %dPs_dFdes{j} = derPs_dFdes(sole,ind_cont3,Ws,dOl_dFdes,dY_dFdes,dFc_dFdes,b);
        %%% dPs_dFdes and dPs_dZdes are the same
        dPs_dZdes{j} = derPs_dZdes(sole,ind_cont3,Ws,dOl_dZdes,dY_dZdes,dFc_dZdes,b);
        dW_Fc_Zdes = Wc * dFc_dZdes;
        dPc_dZdes =  dPs_dZdes{j}(ind_cont3,:);
%         for i=1:3
%             [der_Kctr_dFdes{j}(i,1),der_Kcrot_dFdes{j}(i,1)] = cost_dFdes(dPc_dFdes(:,i),dFc_dFdes(:,i),ind_cont,ind_slip,dOl_dFdes(:,i),dY_dFdes(:,i),ddeltat_dFdes(:,i),Cc,P_mc,Wc,Fc_mc,Fc_mc3,PabsOld_mc,fric,contAngle,angleact,displ_first,Rini,R,Pfree_c,dW_Fc_Fdes(:,i),Xi,Kcart,u_min,u_max,invABwXi,G_inv_A,G_inv_ABw,Fc_mc_hat);   
%         end
        for i=1:2
            [der_Kctr_dZdes{j}(i,1),der_Kcrot_dZdes{j}(i,1)] = cost_dZdes(i,dPc_dZdes(:,i),dFc_dZdes(:,i),ind_cont,ind_slip,dOl_dZdes(:,i),dY_dZdes(:,i),ddeltat_dZdes(:,i),Cc,P_mc,Wc,Fc_mc,Fc_mc3,PabsOld_mc,fric,contAngle,angleact,displ_first,Rini,R,Pfree_c,dW_Fc_Zdes(:,i),Xi,Kcart,u_min,u_max,invABwXi,G_inv_A,G_inv_ABw,Fc_mc_hat,j);   
        end         
    end    
    
    
%     der_Kctr_dZdes = zeros(2,1);
%     der_Kcrot_dZdes = zeros(2,1);
%     derKc_dZdes = cell(2,1);
%     %%% cost_dFdes and cost_dZdes are the same
%     for i=1:2
%         [der_Kctr_dZdes(i),der_Kcrot_dZdes(i),derKc_dZdes{i}] = cost_dZdes(i,dPc_dZdes(:,i),dFc_dZdes(:,i),ind_cont,ind_slip,dOl_dZdes(:,i),dY_dZdes(:,i),ddeltat_dZdes(:,i),Cc,P_mc,Wc,Fc_mc,Fc_mc3,PabsOld_mc,fric,contAngle,angleact,displ_first,Rini,R,Pfree_c,dW_Fc_Zdes(:,i),Xi,Kcart,u_min,u_max,invABwXi,G_inv_A,G_inv_ABw,Fc_mc_hat);   
%     end    
end