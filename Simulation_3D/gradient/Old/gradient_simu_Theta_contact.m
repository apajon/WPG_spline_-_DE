function [Theta_mat,der_HCcH] = gradient_simu_Theta_contact(ind_cont,der_Cc_dp,ind_slip,Fc_mc3,R,Rini,lp,dPFreec_dp,contAngle,dPc_dp)
    der_W_c = zeros(3*length(ind_cont),lp);
%     Hc = zeros(3*length(ind_cont));
%     for j=1:length(ind_cont)
%        r = 3*j-2:3*j;
%        Hc(r,r) = R;
%     end
    dPc_dp(3:3:end,:) = 0;
    der_HCcH = cell(lp,1);
    nc = length(ind_cont);
    if contAngle==1
        t = [1 0 0; 0 1 0];
        Rini1 = [t * Rini;0 0 0];
%         Hc_ini = zeros(3*length(ind_cont));
%         for j=1:length(ind_cont)
%            r = 3*j-2:3*j;
%            Hc_ini(r,r) = Rini1;
%         end    
        for i=1:lp
            tmp = reshape(R*reshape(der_Cc_dp{i},3,3*nc*nc),3*nc,3*nc);
            tmp2 = tmp';
            der_HCcH{i} = reshape(R*reshape(tmp2,3,3*nc*nc),3*nc,3*nc);
            %der_HCcH{i} = Hc * der_Cc_dp{i} * Hc';
            der_HCcHF = der_HCcH{i} * Fc_mc3;            
            der_W_c(:,i) = der_HCcHF + reshape((R - Rini1) * reshape(dPFreec_dp(:,i),3,nc),3*nc,1);
            %der_W_c1 = Hc * der_Cc_dp{i} * Hc' * Fc_mc3 + Hc * dPFreec_dp(:,i) - Hc_ini * dPFreec_dp(:,i);
        end
    else
        for i=1:lp
            tmp = reshape(R*reshape(der_Cc_dp{i},3,3*nc*nc),3*nc,3*nc);
            tmp2 = tmp';
            der_HCcH{i} = reshape(R*reshape(tmp2,3,3*nc*nc),3*nc,3*nc);
            %der_HCcH{i} = Hc * der_Cc_dp{i} * Hc';
            der_HCcHF = der_HCcH{i} * Fc_mc3;          
            %der_W_c(:,i) = der_HCcHF + Hc * dPFreec_dp(:,i) - dPc_dp(:,i);
            der_W_c(:,i) = der_HCcHF + reshape(R* reshape(dPFreec_dp(:,i),3,nc),3*nc,1) - dPc_dp(:,i);
        end
    end
    Theta_mat = [der_W_c;zeros(2*length(ind_slip),lp)];
end