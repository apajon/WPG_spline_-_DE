function dPs_dPfree = derPs_dPfree(sole,angleact,Cs,der_C_s_dPfree,ind_cont,ind_cont3,Ws,Fs3,Pfree_s3,dOl_dPfree,dY_dPfree,dFc_dPfree)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   dPc_dPfree               %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    theta = angleact(1);
    phi = angleact(2);
    psi = angleact(3);
    R = Rot(theta,phi,psi);    
    dPc_dPfree = zeros(3*length(ind_cont),3*sole.nTot);
    dFs_dPfree = zeros(3*sole.nFreeSurf,3*sole.nTot);
    dFs_dPfree(ind_cont3,:) = dFc_dPfree;
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
        for j=1:sole.nFreeSurf
            btheta = dR_dtheta*Pfree_s3((3*j-2):3*j);
            bphi = dR_dphi*Pfree_s3((3*j-2):3*j);
            bpsi = dR_dpsi*Pfree_s3((3*j-2):3*j);
            W_Fs = zeros(3,3);
            RdCRn = zeros(3,1);
            RdCRt1 = zeros(3,1);
            RdCRt2 = zeros(3,1);
            for k=1:sole.nFreeSurf
                btheta = btheta + (dR_dtheta * Cs((3*j-2):3*j,(3*k-2):3*k) * R' + R * Cs((3*j-2):3*j,(3*k-2):3*k) * dR_dtheta') * Fs3((3*k-2):3*k);
                bphi = bphi + (dR_dphi * Cs((3*j-2):3*j,(3*k-2):3*k) * R' + R * Cs((3*j-2):3*j,(3*k-2):3*k) * dR_dphi') * Fs3((3*k-2):3*k);
                bpsi = bpsi + (dR_dpsi * Cs((3*j-2):3*j,(3*k-2):3*k) * R' + R * Cs((3*j-2):3*j,(3*k-2):3*k) * dR_dpsi') * Fs3((3*k-2):3*k);
                W_Fs = W_Fs + Ws((3*j-2):3*j,(3*k-2):3*k) * dFs_dPfree((3*k-2):3*k,(3*g-2):3*g);
                RdCRt1 = RdCRt1 + R * der_C_s_dPfree{(3*g)-2}((3*j-2):3*j,(3*k-2):3*k) * R' * Fs3((3*k-2):3*k);
                RdCRt2 = RdCRt2 + R * der_C_s_dPfree{(3*g)-1}((3*j-2):3*j,(3*k-2):3*k) * R' * Fs3((3*k-2):3*k);
                RdCRn = RdCRn + R * der_C_s_dPfree{3*g}((3*j-2):3*j,(3*k-2):3*k) * R' * Fs3((3*k-2):3*k);
            end
            A = zeros(3,3);
            if ~isempty(find(g==sole.nodesFreeSurf,1))
                if (j==contj && flag==0)
                   A = R;
                   contj = contj + 1;
                   flag = 1;
                end
            end
            dPs_dPfree((3*j-2):3*j,(3*g-2):3*g) = dOl_dPfree(:,(3*g-2):3*g) + [btheta bphi bpsi] * dY_dPfree(:,(3*g-2):3*g) + A + W_Fs + [RdCRt1 RdCRt2 RdCRn];
        end
        flag=0;
        
    end
end