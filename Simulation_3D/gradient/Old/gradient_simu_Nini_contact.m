function Nini = gradient_simu_Nini_contact(ind_cont,Rini,P_mc3,Fc_mc3,angleact,psi_first,contAngle)
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
        PfreeFt1 = 0;
        PfreeFt2 = 0;
        PfreeFn = 0;
        AfreeFt1 = 0;
        AfreeFt2 = 0;
        AfreeFn = 0;
        C = zeros(3,1);
        D = zeros(3,1);
        E = zeros(3,1);
        F = zeros(3,1);
        G = zeros(3,1);
        H = zeros(3,1);
        J = zeros(3,1);
        t1 = [1 0 0];
        t2 = [0 1 0];
        t = [1 0 0; 0 1 0];
        new1 = zeros(1,3);
        new2 = zeros(1,3);
        new3 = zeros(1,3);
        for j=1:length(ind_cont)
            PfreeFt1 = PfreeFt1 + P_mc3((3*j-2):3*j)*Fc_mc3((3*j)-2);
            PfreeFt2 = PfreeFt2 + P_mc3((3*j-2):3*j)*Fc_mc3((3*j)-1);
            PfreeFn = PfreeFn + P_mc3((3*j-2):3*j)*Fc_mc3((3*j));    
%             C = C + (dRinidtheta) * P_mc3((3*j)-2:3*j) * Fc_mc3(3*j);
%             D = D + (dRinidphi) * P_mc3((3*j)-2:3*j) * Fc_mc3(3*j);
%             E = E + (dRinidtheta) * P_mc3((3*j)-2:3*j) * Fc_mc3((3*j)-1);
%             F = F + (dRinidphi) * P_mc3((3*j)-2:3*j) * Fc_mc3((3*j)-1);
%             H = H + (dRinidtheta) * P_mc3((3*j)-2:3*j) * Fc_mc3((3*j)-2);
%             J = J + (dRinidphi) * P_mc3((3*j)-2:3*j) * Fc_mc3((3*j)-2);

% %            
%             n1 = n1 + ([dRinidtheta(2,:)*P_mc3((3*j-2):3*j) dRinidphi(2,:)*P_mc3((3*j-2):3*j) 0]) * Fc_mc3(3*j);
%             n2 = n2 + ([dRinidtheta(1,:)*P_mc3((3*j-2):3*j) dRinidphi(1,:)*P_mc3((3*j-2):3*j) 0]) * Fc_mc3(3*j);
%             n3 = n3 + ([dRinidtheta(1,:)*P_mc3((3*j-2):3*j) dRinidphi(1,:)*P_mc3((3*j-2):3*j) 0]) * Fc_mc3((3*j)-1);
%             n4 = n4 + ([dRinidtheta(2,:)*P_mc3((3*j-2):3*j) dRinidphi(2,:)*P_mc3((3*j-2):3*j) 0]) * Fc_mc3((3*j)-2);
            new1 = new1 - t2 * [dRinidtheta*P_mc3((3*j-2):3*j) dRinidphi*P_mc3((3*j-2):3*j) zeros(3,1)] * Fc_mc3(3*j);
            new2 = new2 + t1 * [dRinidtheta*P_mc3((3*j-2):3*j) dRinidphi*P_mc3((3*j-2):3*j) zeros(3,1)] * Fc_mc3(3*j);
            new3 = new3 - t1 * [dRinidtheta*P_mc3((3*j-2):3*j) dRinidphi*P_mc3((3*j-2):3*j) zeros(3,1)] * Fc_mc3((3*j)-1) + ...
            t2 * [dRinidtheta*P_mc3((3*j-2):3*j) dRinidphi*P_mc3((3*j-2):3*j) zeros(3,1)] * Fc_mc3((3*j)-2);            
% %             PfreeFt1 = PfreeFt1 + P_mc3((3*i-2):3*i)*Fc_mc3((3*i)-2);
%             PfreeFt2 = PfreeFt2 + P_mc3((3*i-2):3*i)*Fc_mc3((3*i)-1);
%             PfreeFn = PfreeFn + P_mc3((3*i-2):3*i)*Fc_mc3((3*i));
        end
        tA = t2 * Rini * Fc_mc3(3*i);
        tB = t1 * Rini * Fc_mc3(3*i);
        tC = -t1 * Rini * Fc_mc3((3*i)-1) + t2 * Rini * Fc_mc3((3*i)-1);
        %A = Rini * Fc_mc3((3*i)-2); 
        %n1 = [-(t2 * C) -(t2 * D) 0];
%         new = zeros(3,3);
%         new(1,1)=dot(dRinidtheta(1,:),PfreeFn);
%         new(1,2)=dot(dRinidphi(1,:),PfreeFn);
%         new(2,1)=-dot(dRinidtheta(2,:),PfreeFn);
%         new(2,2)=-dot(dRinidphi(2,:),PfreeFn);    
%         new(3,1)=-dot(dRinidtheta(1,:),PfreeFt2) + dot(dRinidtheta(2,:),PfreeFt1);
%         new(3,2)=-dot(dRinidphi(1,:),PfreeFt2) + dot(dRinidtheta(2,:),PfreeFt1);
        Nlow(:,(3*i-2):3*i) = [new1; new2; new3] + [tA; tB; tC];

        %Nlow(:,(3*i-2):3*i) = [-(t2 * C) -(t2 * D) 0; (t1 * C) (t1 * D) 0; -(t1 * E)+(t2 * H) -(t1 * F)+(t2 * J) 0] + [tA; tB; tC];
        %Nlow(:,(3*i-2):3*i) = [-n1; n2; (-n3 + n4)];
        %Nlow(:,(3*i-2):3*i) = [-n1' n2' (-n3' + n4')];
    end
    Nini = [zeros(3,3*length(ind_cont));Nlow];    
end

%     Nini = zeros(6,6);
%     theta = angleact(1);
%     phi = angleact(2);
%     psi = angleact(3);    
%     [R,Rtheta,Rphi,Rpsi] = Rot(theta,phi,psi);        
%     Rpsiini = [cos(psi_first), -sin(psi_first), 0; sin(psi_first), cos(psi_first), 0; 0, 0, 1];            
%     dRinidtheta = [0 0 0; 0 -sin(theta) -cos(theta); 0 cos(theta) -sin(theta)]*Rphi*Rpsiini;        
%     dRinidphi = Rtheta*[-sin(phi) 0 cos(phi); 0 0 0; -cos(phi) 0 -sin(phi)]*Rpsiini;    
%     PfreeFt1 = 0;
%     PfreeFt2 = 0;
%     PfreeFn = 0;
%     for i=1:length(ind_cont)
%         PfreeFt1 = PfreeFt1 + P_mc3((3*i-2):3*i)*Fc_mc3((3*i)-2);
%         PfreeFt2 = PfreeFt2 + P_mc3((3*i-2):3*i)*Fc_mc3((3*i)-1);
%         PfreeFn = PfreeFn + P_mc3((3*i-2):3*i)*Fc_mc3((3*i));
%     end
%     Nini(4,4)=dot(Rini(1,:)+dRinidtheta(1,:),PfreeFn);
%     Nini(4,5)=dot(Rini(1,:)+dRinidphi(1,:),PfreeFn);
%     Nini(5,4)=-dot(Rini(2,:)+dRinidtheta(2,:),PfreeFn);
%     Nini(5,5)=-dot(Rini(1,:)+dRinidphi(2,:),PfreeFn);    
%     Nini(6,4)=-dot(Rini(1,:)+dRinidtheta(1,:),PfreeFt2) + dot(dRinidtheta(2,:),PfreeFt1);
%     Nini(6,5)=-dot(Rini(1,:)+dRinidphi(1,:),PfreeFt2) + dot(dRinidtheta(2,:),PfreeFt1);  