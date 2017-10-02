function [dPs_dp,RdPFrees_dp,b] = derPs_dp(sole,angleact,ind_cont,ind_cont3,Ws,Pfree_s3,dOl_dp,dY_dp,dFc_dp,lp,dPFrees_dp,Cs,der_Cs_dp,Fc3)
    theta = angleact(1);
    phi = angleact(2);
    psi = angleact(3);
    R = Rot(theta,phi,psi);    
    dPs_dp = zeros(3*sole.nFreeSurf,lp);
    [dR_dtheta,dR_dphi,dR_dpsi] = dRot(angleact);
    
%     b = cell(sole.nFreeSurf,1);
%     tic;
%     RtF = R'  * reshape(Fc3,3,length(ind_cont));
%     dRthetatF = dR_dtheta' * reshape(Fc3,3,length(ind_cont));
%     dRphitF = dR_dphi' * reshape(Fc3,3,length(ind_cont));
%     dRpsitF = dR_dpsi' * reshape(Fc3,3,length(ind_cont));
%     for j=1:sole.nFreeSurf
%         P = Pfree_s3((3*j-2):3*j);
%         btheta = dR_dtheta*P;
%         bphi = dR_dphi*P;
%         bpsi = dR_dpsi*P;
%         for k=1:length(ind_cont)
%             C = Cs(3*j-2:3*j,3*ind_cont(k)-2:3*ind_cont(k));
%             CRtF = C*RtF(:,k);
%             btheta = btheta + dR_dtheta * CRtF + R * (C * dRthetatF(:,k));
%             bphi = bphi + dR_dphi * CRtF + R * (C * dRphitF(:,k));
%             bpsi = bpsi + dR_dpsi  * CRtF + R * (C * dRpsitF(:,k));
%         end
%         b{j} = [btheta bphi bpsi];
%     end
%     for i=1:lp
%         dPs_dpi = zeros(3*sole.nFreeSurf,1);
%         for j=1:sole.nFreeSurf
%             %P = Pfree_s3((3*j-2):3*j);
%             %btheta = dR_dtheta*P;
%             %bphi = dR_dphi*P;
%             %bpsi = dR_dpsi*P;
%             W_Fc = Ws(3*j-2:3*j,ind_cont3)*dFc_dp(:,i);
%             RdCRt = R * (der_Cs_dp{i}(3*j-2:3*j,ind_cont3)*reshape(RtF,3*length(ind_cont),1));
%             %for k=1:length(ind_cont)
%                 %F = Fc3(3*k-2:3*k);
%                 %C = Cs(3*j-2:3*j,3*ind_cont(k)-2:3*ind_cont(k));
%                 %RtF = R'*F;
%                 %CRtF = C*RtF;
%                 %btheta = btheta + dR_dtheta * CRtF + R * (C * (dR_dtheta' * F));
%                 %bphi = bphi + dR_dphi * CRtF + R * (C * (dR_dphi' * F));
%                 %bpsi = bpsi + dR_dpsi  * CRtF + R * (C * (dR_dpsi' * F));
%                 %W_Fc = W_Fc + Ws(3*j-2:3*j,3*ind_cont(k)-2:3*ind_cont(k)) * dFc_dp((3*k-2):3*k,i);
%                 %RdCRt = RdCRt + R * (der_Cs_dp{i}(3*j-2:3*j,3*ind_cont(k)-2:3*ind_cont(k)) * RtF(:,k));                  
%             %end
%             dPs_dpi(3*j-2:3*j) = dOl_dp(:,i) + b{j} * dY_dp(:,i) + R * dPFrees_dp((3*j-2):3*j,i) + W_Fc + RdCRt;
%         end
%         dPs_dp(:,i) = dPs_dpi;
%     end
%     dPs_dp2=dPs_dp;
    b = zeros(3*sole.nFreeSurf,3);
    %tic;
    RtF = R'  * reshape(Fc3,3,length(ind_cont));
    dRthetatF = dR_dtheta' * reshape(Fc3,3,length(ind_cont));
    dRphitF = dR_dphi' * reshape(Fc3,3,length(ind_cont));
    dRpsitF = dR_dpsi' * reshape(Fc3,3,length(ind_cont));
    for j=1:sole.nFreeSurf
        P = Pfree_s3((3*j-2):3*j);
%         btheta = dR_dtheta*P;
%         bphi = dR_dphi*P;
%         bpsi = dR_dpsi*P;
        btheta = dR_dtheta*P + R * (Cs(3*j-2:3*j,ind_cont3)*reshape(dRthetatF,3*length(ind_cont),1))...
                         + dR_dtheta * (Cs(3*j-2:3*j,ind_cont3)*reshape(RtF,3*length(ind_cont),1));
        bphi = dR_dphi*P + R * (Cs(3*j-2:3*j,ind_cont3)*reshape(dRphitF,3*length(ind_cont),1))...
                     + dR_dphi * (Cs(3*j-2:3*j,ind_cont3)*reshape(RtF,3*length(ind_cont),1));
        bpsi = dR_dpsi*P + R * (Cs(3*j-2:3*j,ind_cont3)*reshape(dRpsitF,3*length(ind_cont),1))...
                     + dR_dpsi * (Cs(3*j-2:3*j,ind_cont3)*reshape(RtF,3*length(ind_cont),1));
%         for k=1:length(ind_cont)
%             C = Cs(3*j-2:3*j,3*ind_cont(k)-2:3*ind_cont(k));
%             CRtF = C*RtF(:,k);
%             btheta = btheta + dR_dtheta * CRtF + R * (C * dRthetatF(:,k));
%             bphi = bphi + dR_dphi * CRtF + R * (C * dRphitF(:,k));
%             bpsi = bpsi + dR_dpsi  * CRtF + R * (C * dRpsitF(:,k));
%         end
        b(3*j-2:3*j,:) = [btheta bphi bpsi];
    end
    
    W_Fc = Ws(:,ind_cont3)*dFc_dp;
    for i=1:lp
        %dPs_dpi = zeros(3*sole.nFreeSurf,1);
        %for j=1:sole.nFreeSurf
            %P = Pfree_s3((3*j-2):3*j);
            %btheta = dR_dtheta*P;
            %bphi = dR_dphi*P;
            %bpsi = dR_dpsi*P;
            
            
            
            
            %RdCRt = R * (der_Cs_dp{i}(3*j-2:3*j,ind_cont3)*reshape(RtF,3*length(ind_cont),1));
            
        tmp = der_Cs_dp{i}(:,ind_cont3)*reshape(RtF,3*length(ind_cont),1);
        RdCRt = reshape(R*reshape(tmp,3,sole.nFreeSurf),3*sole.nFreeSurf,1);
            %for k=1:length(ind_cont)
                %F = Fc3(3*k-2:3*k);
                %C = Cs(3*j-2:3*j,3*ind_cont(k)-2:3*ind_cont(k));
                %RtF = R'*F;
                %CRtF = C*RtF;
                %btheta = btheta + dR_dtheta * CRtF + R * (C * (dR_dtheta' * F));
                %bphi = bphi + dR_dphi * CRtF + R * (C * (dR_dphi' * F));
                %bpsi = bpsi + dR_dpsi  * CRtF + R * (C * (dR_dpsi' * F));
                %W_Fc = W_Fc + Ws(3*j-2:3*j,3*ind_cont(k)-2:3*ind_cont(k)) * dFc_dp((3*k-2):3*k,i);
                %RdCRt = RdCRt + R * (der_Cs_dp{i}(3*j-2:3*j,3*ind_cont(k)-2:3*ind_cont(k)) * RtF(:,k));                  
            %end
            
            %dPs_dpi(3*j-2:3*j) = 
        %end
        RdPFrees_dp = reshape(R*reshape(dPFrees_dp(:,i),3,sole.nFreeSurf),3*sole.nFreeSurf,1);
        dPs_dp(:,i) = repmat(dOl_dp(:,i),sole.nFreeSurf,1) + b * dY_dp(:,i) + RdPFrees_dp + W_Fc(:,i) + RdCRt;
    end
    %toc
end