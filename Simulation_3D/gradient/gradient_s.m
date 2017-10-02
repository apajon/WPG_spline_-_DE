function [G,A,Bw,Xi,Ws,Hc,inv_A,E,A_1B,inv_E_A1B,u_max,u_min] = gradient_s(sole,Kcart,displ,angleact,Zdes,fric,Pfree_s,PSurf3,PabsOld,Fc_mc3,Cs,ind_cont,ind_slip,displ_first,psi_first,s)
    ind_cont3 = sort(([3*ind_cont-2; 3*ind_cont-1; 3*ind_cont]));
    ind_cont3 = reshape(ind_cont3,3*length(ind_cont),1);  
    theta = angleact(1);
    phi = angleact(2);
    psi = angleact(3);
    R = Rot(theta,phi,psi);
    Rini = Rot(theta,phi,psi_first);
    Pfree_c = Pfree_s(:,ind_cont);
    Pfree_c3 = reshape(Pfree_c,3*length(ind_cont),1);
    %Pfree_s3 = reshape(Pfree_s,3*sole.nFreeSurf,1);
    Fc_mc = reshape(Fc_mc3,3,length(ind_cont));
    P_mc3 = PSurf3(ind_cont3);
    P_mc = reshape(P_mc3,3,length(ind_cont)); 
    Ws = zeros(3*sole.nFreeSurf,3*sole.nFreeSurf);
    for j = 1:sole.nFreeSurf
        for k = 1:sole.nFreeSurf
            Ws(3*j-2:3*j,(3*k)-2:3*k) = R * Cs(3*j-2:3*j,3*k-2:3*k) * R';
        end
    end
    Cc = Cs(ind_cont3,ind_cont3);
    Wc = Ws(ind_cont3,ind_cont3);
    PabsOld_mc = PabsOld(:,ind_cont);
    %dPFrees_dp = dPFree_dp(sole.nodesFreeSurf3,:);
    %dPFreec_dp = dPFrees_dp(ind_cont3,:);
    %%% Compute Kcart     
    A = gradient_simu_A_contact(P_mc,Fc_mc,ind_cont,Wc,ind_slip,PabsOld_mc,fric,s,Rini,displ_first);
    G = computeG(ind_cont,ind_slip,Fc_mc,P_mc,Zdes);
    Bw = computeBw(Fc_mc3,ind_cont,ind_slip,R,Wc,Pfree_c);
    Xi = computeXi(displ,Zdes);
    %Kcart = computeKc(G,A,Bw,Xi);    
    %%% Compute gradient of Kcart
%     if ~isempty(dPs_dp_previous)
%         dPc_dp_previous = dPs_dp_previous(ind_cont3,:);
%     else
%         dPc_dp_previous = [];
%     end
    [B,B12out] = gradient_simu_B_contact(angleact,Pfree_c3,Fc_mc3,ind_cont,Cc,ind_slip,psi_first,s);
    E = gradient_simu_E_contact(ind_cont,ind_slip,Fc_mc,P_mc,Zdes);
    inv_A = inv(A);
    A_1B = inv_A * B;
    E_A1B = E * A_1B;

    if s==1
       N_ini = gradient_simu_N_ini_contact(ind_cont,Pfree_c3,Fc_mc3,angleact,psi_first,s);
       inv_E_A1B = inv(E_A1B + N_ini);
    else
       inv_E_A1B = inv(E_A1B);
    end
    Hc = zeros(3*length(ind_cont));
    for j=1:length(ind_cont)
        r = 3*j-2:3*j;
        Hc(r,r) = R;
    end     
    Kc_squared_tr = Kcart(1:3,1:3)'*Kcart(1:3,1:3);
    [u,DE] = eig(Kc_squared_tr);
    [~,I_lamba_max] = max(diag(DE));
    u_max = u(:,I_lamba_max);
    Kc_squared_rot = Kcart(4:6,4:6)'*Kcart(4:6,4:6);
    [u,DE] = eig(Kc_squared_rot);
    [~,I_lamba_min] = min(diag(DE));
    u_min = u(:,I_lamba_min);
end