function [Theta_mat_pi,der_HCcH_pi,der_HCcHF_pi] = gradient_simu_Theta_contact_pi(ind_cont,der_Cc_dpi,ind_slip,Fc_mc3,Hc,Rini,dPFreec_dpi,contAngle,dPc_dpi_previous)  
    if contAngle==1
        t = [1 0 0; 0 1 0];
        Rini1 = [t * Rini;0 0 0];
        Hc_ini = zeros(3*length(ind_cont));
        for j=1:length(ind_cont)
            r = 3*j-2:3*j;
            Hc_ini(r,r) = Rini1;
        end    
        der_HCcH_pi = Hc * der_Cc_dpi * Hc';
        der_HCcHF_pi = der_HCcH_pi * Fc_mc3;            
        Theta_mat_pi = [der_HCcHF_pi + (Hc - Hc_ini) * dPFreec_dpi;zeros(2*length(ind_slip),1)];
    else
        dPc_dpi_previous(3:3:end,:) = 0;
        der_HCcH_pi = Hc * der_Cc_dpi * Hc';
        der_HCcHF_pi = der_HCcH_pi * Fc_mc3;          
        Theta_mat_pi = [der_HCcHF_pi + Hc * dPFreec_dpi - dPc_dpi_previous;zeros(2*length(ind_slip),1)];
    end
end