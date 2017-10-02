function [dPs_dp,dW_Fc] = derPs_dp(sole,angleact,ind_cont,Wc,Pfree_s3,dOl_dp,dY_dp,dFc_dp,lp,dPFrees_dp,B12out,der_HCcHF)
    theta = angleact(1);
    phi = angleact(2);
    psi = angleact(3);
    R = Rot(theta,phi,psi);    
    dPs_dp = zeros(3*sole.nFreeSurf,lp);
    [dR_dtheta,dR_dphi,dR_dpsi] = dRot(angleact);
    dW_Fc = Wc * dFc_dp;
    
    cont = 0;
    for j=1:sole.nFreeSurf
        if any(j==ind_cont)
            cont = cont + 1;
            b = B12out((3*cont-2):3*cont,:);
            dPs_dp((3*j-2):3*j,:) = dOl_dp + b * dY_dp + R * dPFrees_dp((3*j-2):3*j,:) + dW_Fc(3*cont-2:3*cont,:);
            for g=1:lp
                 dPs_dp((3*j-2):3*j,g)  = dPs_dp((3*j-2):3*j,g) + der_HCcHF{g}(3*cont-2:3*cont);
            end
            %for k=1:length(ind_cont)
                %F = Fc3((3*k-2):3*k);
                %C = Cc((3*cont-2):3*cont,(3*k-2):3*k);
                %RtF = R'*F;
%                     CRtF = C*RtF;
%                     btheta = btheta + dR_dtheta * CRtF + R * (C * (dR_dtheta' * F));
%                     bphi = bphi + dR_dphi * CRtF + R * (C * (dR_dphi' * F));
%                     bpsi = bpsi + dR_dpsi  * CRtF + R * (C * (dR_dpsi' * F));
                %W_Fc = W_Fc + Wc((3*cont-2):3*cont,(3*k-2):3*k) * dFc_dp((3*k-2):3*k,g);
                %RdCRt = RdCRt + R * (der_Cc_dp{g}((3*cont-2):3*cont,(3*k-2):3*k) * RtF);                  
            %end
        else 
            P = Pfree_s3((3*j-2):3*j);
            btheta = dR_dtheta*P;
            bphi = dR_dphi*P;
            bpsi = dR_dpsi*P;
            b = [btheta bphi bpsi];
            dPs_dp((3*j-2):3*j,:) = dOl_dp + b * dY_dp + R * dPFrees_dp((3*j-2):3*j,:);
        end
    end
end