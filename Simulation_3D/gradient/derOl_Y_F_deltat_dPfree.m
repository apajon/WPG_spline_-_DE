function [dOl_dPfree,dY_dPfree,dFc_dPfree,ddeltat_dPfree] = derOl_Y_F_deltat_dPfree(sole,Zdes,angleact,Ccc,der_C_cc_dPfree,ind_cont,ind_slip,Wc,Fc_mc,Fc_mc3,P_mc,Pfree_c,Pfree_c3,PabsOld_mc,fric,contAngle,displ_first,psi_first)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  Ol_dPfree, dY_dPfree and dF_dPfree%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    theta = angleact(1);
    phi = angleact(2);
    psi = angleact(3);
    R = Rot(theta,phi,psi);
    Rini = Rot(theta,phi,psi_first);
    A = gradient_simu_A_contact(P_mc,Fc_mc,ind_cont,Wc,ind_slip,Pfree_c,fric,contAngle,Rini,displ_first);
    B = gradient_simu_B_contact(angleact,Pfree_c3,Fc_mc3,ind_cont,Ccc,ind_slip,psi_first,contAngle);
    E = gradient_simu_E_contact(ind_cont,ind_slip,Fc_mc,P_mc,Zdes);
    N_ini = gradient_simu_N_ini_contact(ind_cont,Pfree_c3,Fc_mc3,angleact,psi_first,contAngle);
    Theta_mat = gradient_simu_Theta_contact(sole,ind_cont,der_C_cc_dPfree,ind_slip,Fc_mc3,R,Rini);
    inv_A = inv(A);
    A_1B = inv_A * B;
    E_A1B = (E * A_1B);
    E_AN_ini = inv(E_A1B + N_ini);
    
    inv_A_theta = inv_A * Theta_mat;
    EinvTheta = -E * inv_A_theta;
    dOl_dY = E_AN_ini * EinvTheta;
    contact_surf = 0;
    for j=1:sole.nTot
        if ~isempty(find(j==sole.nodesFreeSurf,1))
            contact_surf = contact_surf + 1;
            contact = find(contact_surf==ind_cont,1);
        else
            contact = 0;
        end
        if ~isempty(find(contact_surf==ind_cont,1))
            Y = gradient_simu_Y_contact(contact,Rini,Fc_mc3);             
            dOl_dY(:,(3*j)-2:(3*j)) = dOl_dY(:,(3*j)-2:(3*j)) + E_AN_ini * Y;
%             dF_ddelta_dPfree = A_1B * dOl_dY(:,(3*j)-2:(3*j)) + inv_A_theta;
%         else
%             dOl_dY(:,(3*j)-2:(3*j)) = E_AN_ini * EinvTheta;
%             dF_ddelta_dPfree = A_1B * dOl_dY(:,(3*j)-2:(3*j)) + inv_A_theta;
%             dFc_dPfree(:,(3*j)-2:(3*j)) = dF_ddelta_dPfree(1:3*length(ind_cont),:);    
%             ddeltat_dPfree(:,(3*j)-2:(3*j)) = dF_ddelta_dPfree(3*length(ind_cont)+1:end,:);
         end
    end
    dF_ddelta_dPfree = A_1B * dOl_dY + inv_A_theta;
    dFc_dPfree = dF_ddelta_dPfree(1:3*length(ind_cont),:);
    ddeltat_dPfree = dF_ddelta_dPfree(3*length(ind_cont)+1:end,:);

    dOl_dPfree = dOl_dY(1:3,:);
    dY_dPfree = dOl_dY(4:6,:); 
end