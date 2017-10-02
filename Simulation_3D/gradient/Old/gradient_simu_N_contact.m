function N = gradient_simu_N_contact(ind_cont,Rini,P_mc3,Fc_mc3,angleact,psi_first,contAngle)
    [dR_dtheta,dR_dphi,dR_dpsi] = dRot(angleact);
    t2 = [0 1 0];
    t1 = [1 0 0];
    if contAngle==1
        theta = angleact(1);
        phi = angleact(2);
        psi = angleact(3);    
        [R,Rtheta,Rphi,Rpsi] = Rot(theta,phi,psi);        
        Rpsiini = [cos(psi_first), -sin(psi_first), 0; sin(psi_first), cos(psi_first), 0; 0, 0, 1];            
        dRinidtheta = [0 0 0; 0 -sin(theta) -cos(theta); 0 cos(theta) -sin(theta)]*Rphi*Rpsiini;        
        dRinidphi = Rtheta*[-sin(phi) 0 cos(phi); 0 0 0; -cos(phi) 0 -sin(phi)]*Rpsiini;
    end
    Nlow = zeros(3,3*length(ind_cont));
    for i=1:length(ind_cont)
        n1 = zeros(1,3);
        n2 = zeros(1,3);
        n3 = zeros(1,3);
        n4 = zeros(1,3);
        btheta = zeros(3,1);
        bphi = zeros(3,1);
        for j=1:length(ind_cont)
            if i==j
                A = Rini * [1 0 0;0 1 0; 0 0 1];
            else
                A = zeros(3,3);
            end
%             %btheta = -dRinidtheta*P_mc3((3*j-2):3*j);
%             %bphi = -dRinidphi*P_mc3((3*j-2):3*j); 
%             %????????????????????????????
%             %btheta = dRinidtheta*P_mc3((3*j-2):3*j);
%             %bphi = dRinidphi*P_mc3((3*j-2):3*j);             
%             btheta = dRinidtheta*P_mc3((3*j-2):3*j);
%             bphi = dRinidphi*P_mc3((3*j-2):3*j);
%             %B1 = [btheta'; bphi'; zeros(1,3)];
%             B1 = [btheta bphi zeros(3,1)];
%             n1 = n1 + t2 * (A + B1) * Fc_mc3(3*j);
%             n2 = n2 + t1 * (A + B1) * Fc_mc3(3*j);
%             n3 = n3 + t1 * (A + B1) * Fc_mc3((3*j)-1);
%             n4 = n4 + t2 * (A + B1) * Fc_mc3((3*j)-2);
            n1 = n1 + (A(2,:) + [dRinidtheta(2,:)*P_mc3((3*j-2):3*j) dRinidphi(2,:)*P_mc3((3*j-2):3*j) 0]) * Fc_mc3(3*j);
            n2 = n2 + (A(1,:) + [dRinidtheta(1,:)*P_mc3((3*j-2):3*j) dRinidphi(1,:)*P_mc3((3*j-2):3*j) 0]) * Fc_mc3(3*j);
            n3 = n3 + (A(1,:) + [dRinidtheta(1,:)*P_mc3((3*j-2):3*j) dRinidphi(1,:)*P_mc3((3*j-2):3*j) 0]) * Fc_mc3((3*j)-1);
            n4 = n4 + (A(2,:) + [dRinidtheta(2,:)*P_mc3((3*j-2):3*j) dRinidphi(2,:)*P_mc3((3*j-2):3*j) 0]) * Fc_mc3((3*j)-2);
        end
        Nlow(:,(3*i-2):3*i) = [-n1; n2; (-n3 + n4)];
        %Nlow(:,(3*i-2):3*i) = [-n1' n2' (-n3' + n4')];
    end
    N = [zeros(3,3*length(ind_cont));Nlow];   
end