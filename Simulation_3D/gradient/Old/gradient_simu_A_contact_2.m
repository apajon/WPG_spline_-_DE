function A = gradient_simu_A_contact_2(P_mc,Fc_mc,ind_cont,Wc,ind_slip,PabsOld_mc,fric,contAngle,displ_ini,psi_ini,angleact)
    A11 = -Wc;
    A12 = zeros(3*length(ind_cont),2*length(ind_slip));
    A21 = zeros(2*length(ind_slip),3*length(ind_cont));
    A22 = zeros(2*length(ind_slip),2*length(ind_slip));
    if ~isempty(ind_slip)
        Islip = [];
        for i=1:length(ind_slip)
            Islip(i) = find(ind_cont==ind_slip(i));
        end
        Fc_mc_splip = Fc_mc(:,Islip);
        for i=1:length(ind_slip)
            if contAngle==1
                Rpsiini = [cos(psi_ini), -sin(psi_ini), 0; sin(psi_ini), cos(psi_ini), 0; 0, 0, 1];
                theta = angleact(1);
                phi = angleact(2);
                psi = angleact(3);
                [R,Rtheta,Rphi,Rpsi] = Rot(theta,phi,psi);
                Rini = Rtheta * Rphi * Rpsiini;
                delta = P_mc(:,Islip(i)) - displ_ini - Rini * PabsOld_mc(:,Islip(i));
            else
                delta = P_mc(:,Islip(i))-PabsOld_mc(:,Islip(i));
            end
            delta_t = delta(1:2,1);
            A21((2*i)-1:2*i,(3*Islip(i))-2:3*Islip(i)) = [norm(delta_t).*eye(2) fric.*delta_t];
            A12((3*Islip(i))-2:3*Islip(i),(2*i)-1:2*i) = [1 0; 0 1; 0 0];
            A22((2*i)-1:2*i,(2*i)-1:2*i) =  [1;1]*(delta_t'./norm(delta_t));
        end
    end
    A = [A11 A12;A21 A22];
end