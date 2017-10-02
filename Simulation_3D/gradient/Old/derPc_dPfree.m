function [dPc_dPfree] = derPc_dPfree(sole,angleact,Ccc,der_C_cc_dPfree,ind_cont,Wc,Fc_mc3,Pfree_c3,dOl_dPfree,dY_dPfree,dFc_dPfree)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   dPc_dPfree               %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    theta = angleact(1);
    phi = angleact(2);
    psi = angleact(3);
    R = Rot(theta,phi,psi);    
    dPc_dPfree = zeros(3*length(ind_cont),3*sole.nTot);
    [dR_dtheta,dR_dphi,dR_dpsi] = dRot(angleact);
    contj = 1;
    flag=0;
    contact_surf = 0;
    for g=1:sole.nTot
        if ~isempty(find(g==sole.nodesFreeSurf,1))
            contact_surf = contact_surf + 1;
            contact = find(contact_surf==ind_cont,1);
        else
            contact = 0;
        end 
        for j=1:length(ind_cont)
            btheta = dR_dtheta*Pfree_c3((3*j-2):3*j);
            bphi = dR_dphi*Pfree_c3((3*j-2):3*j);
            bpsi = dR_dpsi*Pfree_c3((3*j-2):3*j);
            W_Fs = zeros(3,3);
            RdCRn = zeros(3,1);
            RdCRt1 = zeros(3,1);
            RdCRt2 = zeros(3,1);
            for k=1:length(ind_cont)
                btheta = btheta + (dR_dtheta * Ccc((3*j-2):3*j,(3*k-2):3*k) * R' + R * Ccc((3*j-2):3*j,(3*k-2):3*k) * dR_dtheta') * Fc_mc3((3*k-2):3*k);
                bphi = bphi + (dR_dphi * Ccc((3*j-2):3*j,(3*k-2):3*k) * R' + R * Ccc((3*j-2):3*j,(3*k-2):3*k) * dR_dphi') * Fc_mc3((3*k-2):3*k);
                bpsi = bpsi + (dR_dpsi * Ccc((3*j-2):3*j,(3*k-2):3*k) * R' + R * Ccc((3*j-2):3*j,(3*k-2):3*k) * dR_dpsi') * Fc_mc3((3*k-2):3*k);
                W_Fs = W_Fs + Wc((3*j-2):3*j,(3*k-2):3*k) * dFc_dPfree((3*k-2):3*k,(3*g-2):3*g); 
                RdCRt1 = RdCRt1 + R * der_C_cc_dPfree{(3*g)-2}((3*j-2):3*j,(3*k-2):3*k) * R' * Fc_mc3((3*k-2):3*k);
                RdCRt2 = RdCRt2 + R * der_C_cc_dPfree{(3*g)-1}((3*j-2):3*j,(3*k-2):3*k) * R' * Fc_mc3((3*k-2):3*k);
                RdCRn = RdCRn + R * der_C_cc_dPfree{3*g}((3*j-2):3*j,(3*k-2):3*k) * R' * Fc_mc3((3*k-2):3*k);
            end
            A = zeros(3,3);
            if ~isempty(find(g==sole.nodesFreeSurf,1))
                if ~isempty(find(contact_surf==ind_cont,1))
                    if (j==contj && flag==0)
                       A = R;
                       contj = contj + 1;
                       flag = 1;
                    end
                end
            end
            dPc_dPfree((3*j-2):3*j,(3*g-2):3*g) = dOl_dPfree(:,(3*g-2):3*g) + [btheta bphi bpsi] * dY_dPfree(:,(3*g-2):3*g) + A + W_Fs + [RdCRt1 RdCRt2 RdCRn];
        end
        flag=0;
        
    end
end