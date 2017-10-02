function B = gradient_simu_B_contact(angleact,Pfree_c3,Fc_mc3,ind_cont,Cc,ind_slip,psi_first,contAngle)
    [dR_dtheta,dR_dphi,dR_dpsi] = dRot(angleact);
    theta = angleact(1);
    phi = angleact(2);
    psi = angleact(3);    
    [R,Rtheta,Rphi,Rpsi] = Rot(theta,phi,psi);
    if contAngle==1
        Rpsiini = [cos(psi_first) -sin(psi_first) 0; sin(psi_first) cos(psi_first) 0; 0 0 1];            
        dRinidtheta = [0 0 0; 0 -sin(theta) -cos(theta); 0 cos(theta) -sin(theta)]*Rphi*Rpsiini;        
        dRinidphi = Rtheta*[-sin(phi) 0 cos(phi); 0 0 0; -cos(phi) 0 -sin(phi)]*Rpsiini;
    end
    for i=1:length(ind_cont)
        if contAngle==1
            btheta = dR_dtheta*Pfree_c3((3*i-2):3*i);
            bphi = dR_dphi*Pfree_c3((3*i-2):3*i);
            bpsi = dR_dpsi*Pfree_c3((3*i-2):3*i);
            btheta(1) = btheta(1) - dRinidtheta(1,:)*Pfree_c3((3*i-2):3*i);
            btheta(2) = btheta(2) - dRinidtheta(2,:)*Pfree_c3((3*i-2):3*i);
            bphi(1) = bphi(1) - dRinidphi(1,:)*Pfree_c3((3*i-2):3*i);
            bphi(2) = bphi(2) - dRinidphi(2,:)*Pfree_c3((3*i-2):3*i);
        else
            btheta = dR_dtheta*Pfree_c3((3*i-2):3*i);
            bphi = dR_dphi*Pfree_c3((3*i-2):3*i);
            bpsi = dR_dpsi*Pfree_c3((3*i-2):3*i);           
        end
        for j=1:length(ind_cont)
            F = Fc_mc3(3*j-2:3*j);
            C = Cc(3*i-2:3*i,3*j-2:3*j);
            RtF = R'*F;
            CRtF = C*RtF;
            a = dR_dtheta * CRtF + R * (C * (dR_dtheta' * F));
            b = dR_dphi * CRtF + R * (C * (dR_dphi' * F));
            c = dR_dpsi  * CRtF + R * (C * (dR_dpsi' * F));
            btheta = btheta + a;
            bphi = bphi + b; 
            bpsi = bpsi + c;
%             btheta = btheta + (dR_dtheta * Ccc((3*i-2):3*i,(3*j-2):3*j) * R' + R * Ccc((3*i-2):3*i,(3*j-2):3*j) * dR_dtheta') * Fc_mc3((3*j-2):3*j);
%             bphi = bphi + (dR_dphi * Ccc((3*i-2):3*i,(3*j-2):3*j) * R' + R * Ccc((3*i-2):3*i,(3*j-2):3*j) * dR_dphi') * Fc_mc3((3*j-2):3*j);
%             bpsi = bpsi + (dR_dpsi * Ccc((3*i-2):3*i,(3*j-2):3*j) * R' + R * Ccc((3*i-2):3*i,(3*j-2):3*j) * dR_dpsi') * Fc_mc3((3*j-2):3*j);
        end
        B12((3*i-2):3*i,:) = [btheta bphi bpsi];
    end
    B11 = repmat(eye(3),length(ind_cont),1);
    B = [B11 B12;zeros(2*length(ind_slip),3) zeros(2*length(ind_slip),3)];
end