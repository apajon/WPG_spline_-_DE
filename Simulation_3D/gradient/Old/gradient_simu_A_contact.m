function A = gradient_simu_A_contact(P_mc,Fc_mc,ind_cont,Wc,ind_slip,PabsOld_mc,fric,contAngle,Rini,displ_first)
    t = [1 0 0; 0 1 0];
    A11 = -Wc;
    A12 = zeros(3*length(ind_cont),2*length(ind_slip));
    A21 = zeros(2*length(ind_slip),3*length(ind_cont));
    A22 = zeros(2*length(ind_slip),2*length(ind_slip));
    if ~isempty(ind_slip)
        Islip = zeros(length(ind_slip),1);
        for i=1:length(ind_slip)
            Islip(i) = find(ind_cont==ind_slip(i));
        end
        Fc_mc_splip = Fc_mc(:,Islip);
        P_mc_slip = P_mc(:,Islip);
        PabsOld_mc_slip = PabsOld_mc(:,Islip);
        for i=1:length(ind_slip)
            if contAngle==1
                delta = P_mc_slip(:,i) - [t * displ_first + t * Rini * PabsOld_mc_slip(:,i);0];
            else
                delta = P_mc_slip(:,i) - [t * PabsOld_mc_slip(:,i);0];
            end
            delta_t = delta(1:2,1);
            A21((2*i)-1:2*i,(3*Islip(i))-2:3*Islip(i)) = [norm(delta_t)*eye(2) fric*delta_t];
            A12((3*Islip(i))-2:3*Islip(i),(2*i)-1:2*i) = [1 0; 0 1; 0 0];
            A22((2*i)-1:2*i,(2*i)-1:2*i) =  Fc_mc_splip(1:2,i)*(delta_t'./norm(delta_t))+fric*Fc_mc_splip(3,i)*eye(2);
        end
    end
    A = [A11 A12;A21 A22];
end