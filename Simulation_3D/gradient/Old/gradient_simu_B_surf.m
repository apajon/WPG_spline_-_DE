function B = gradient_simu_B_surf(sole,angleact,Pfree_s3,FSurf3,Cs,ind_slip,psi_first,contAngle)
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
    for i=1:sole.nFreeSurf
        if contAngle==1
            btheta = dR_dtheta*Pfree_s3((3*i-2):3*i);
            bphi = dR_dphi*Pfree_s3((3*i-2):3*i);
            bpsi = dR_dpsi*Pfree_s3((3*i-2):3*i);
            btheta(1) = btheta(1) - dRinidtheta(1,:)*Pfree_s3((3*i-2):3*i);
            btheta(2) = btheta(2) - dRinidtheta(2,:)*Pfree_s3((3*i-2):3*i);
            bphi(1) = bphi(1) - dRinidphi(1,:)*Pfree_s3((3*i-2):3*i);
            bphi(2) = bphi(2) - dRinidphi(2,:)*Pfree_s3((3*i-2):3*i);
        else
            btheta = dR_dtheta*Pfree_s3((3*i-2):3*i);
            bphi = dR_dphi*Pfree_s3((3*i-2):3*i);
            bpsi = dR_dpsi*Pfree_s3((3*i-2):3*i);
        end
        for j=1:sole.nFreeSurf
            btheta = btheta + (dR_dtheta * Cs((3*i-2):3*i,(3*j-2):3*j) * R' + R * Cs((3*i-2):3*i,(3*j-2):3*j) * dR_dtheta') * FSurf3((3*j-2):3*j);
            bphi = bphi + (dR_dphi * Cs((3*i-2):3*i,(3*j-2):3*j) * R' + R * Cs((3*i-2):3*i,(3*j-2):3*j) * dR_dphi') * FSurf3((3*j-2):3*j);
            bpsi = bpsi + (dR_dpsi * Cs((3*i-2):3*i,(3*j-2):3*j) * R' + R * Cs((3*i-2):3*i,(3*j-2):3*j) * dR_dpsi') * FSurf3((3*j-2):3*j);
        end
        B12((3*i-2):3*i,:) = [btheta bphi bpsi];
    end
    B11 = repmat(eye(3),sole.nFreeSurf,1);
    B = [B11 B12;zeros(2*length(ind_slip),3) zeros(2*length(ind_slip),3)];
end