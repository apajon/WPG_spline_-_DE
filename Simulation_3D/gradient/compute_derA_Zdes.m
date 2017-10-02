function [derA_Zdes,derWab_i,dR_dZdes]  = compute_derA_Zdes(dY_dZdes,dFc_dZdes,ddeltat_dZdes,Ccc,P_mc,Fc_mc,ind_cont,ind_slip,PabsOld_mc,fric,contAngle,angleact,displ_first,Rini,R)
    t = [1 0 0; 0 1 0];
    [dR_dtheta,dR_dphi,dR_dpsi] = dRot(angleact);
    derA12_Zdes = zeros(3*length(ind_cont),2*length(ind_slip));
    derA21_Zdes = zeros(2*length(ind_slip),3*length(ind_cont));
    derA22_Zdes = zeros(2*length(ind_slip),2*length(ind_slip));
    dR_dZdes = dR_dtheta .* dY_dZdes(1) + dR_dphi .* dY_dZdes(2) + dR_dpsi .* dY_dZdes(3);
    nc = length(ind_cont);
    tmp = reshape(R*reshape(Ccc,3,3*nc*nc),3*nc,3*nc);
    tmp2 = tmp';
    dR_1 = reshape(dR_dZdes*reshape(tmp2,3,3*nc*nc),3*nc,3*nc);
    derA11_Zdes = dR_1 + dR_1';
    if ~isempty(ind_slip)
        Islip = zeros(length(ind_slip),1);
        for i=1:length(ind_slip)
            Islip(i) = find(ind_cont==ind_slip(i));
        end
        Fc_mc_splip = Fc_mc(:,Islip);
        P_mc_slip = P_mc(:,Islip);
        PabsOld_mc_slip = PabsOld_mc(:,Islip);
        for i=1:length(ind_slip)
            if contAngle==1
                delta = P_mc_slip(:,i) - [t * displ_first + t * Rini * PabsOld_mc_slip(:,i);0];
            else
                delta = P_mc_slip(:,i) - [t * PabsOld_mc_slip(:,i);0];
            end
            delta_t = delta(1:2);
            derA21_Zdes((2*i)-1:2*i,(3*Islip(i))-2:3*Islip(i)) = [(delta_t)'./norm(delta_t)*ddeltat_dZdes(2*i-1:2*i)*eye(2) fric.*ddeltat_dZdes(2*i-1:2*i)];
            a1 = dFc_dZdes(3*Islip(i)-2:3*Islip(i)-1)*((delta_t)'./norm(delta_t));
            a2 = Fc_mc_splip(1:2,i) *(((eye(2).*norm(delta_t)-delta_t*(delta_t'./norm(delta_t)))./(norm(delta_t)^2))*ddeltat_dZdes(2*i-1:2*i))';
            derA22_Zdes((2*i)-1:2*i,(2*i)-1:2*i) = a1 + a2 + fric * dFc_dZdes(3*Islip(i))*eye(2);
        end
    end
    derWab_i = derA11_Zdes;
    derA_Zdes = [-derA11_Zdes derA12_Zdes;derA21_Zdes derA22_Zdes];
end