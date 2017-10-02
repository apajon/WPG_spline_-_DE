function derG_p_i = compute_derG_pi(dPc_dp_i,dFc_dp_i,ind_cont,ind_slip)
%     contact_nTot = [];
    derG21_pi = zeros(3,3*length(ind_cont));
%     cont = 1;
%     contact_surf = 0;
%     for j=1:sole.nTot
%         if ~isempty(find(j==sole.nodesFreeSurf,1))
%             contact_surf = contact_surf + 1;
%             if ~isempty(find(contact_surf==ind_cont,1))
%                 contact_nTot(cont) = j;
%                 cont = cont + 1;
%             end
%         end
%     end
    for i=1:length(ind_cont)
        derG21_pi(:,3*i-2:3*i) = [0 0 dPc_dp_i(3*i-1);0 0 -dPc_dp_i(3*i-2); -dPc_dp_i(3*i-1) dPc_dp_i(3*i-2) 0]; 
    end
    derG22_pi = [];
    if ~isempty(ind_slip)      
        Islip = zeros(length(ind_slip),1);
        for i=1:length(ind_slip)
            Islip(i) = find(ind_cont==ind_slip(i));
        end
        derG22_pi = [];
        for i=1:length(ind_slip)
            derG22_pi(:,2*i-1:2*i) = [0 dFc_dp_i(3*Islip(i));dFc_dp_i(3*Islip(i)) 0; dFc_dp_i(3*Islip(i)-1) -dFc_dp_i(3*Islip(i)-2)]; 
        end    
    end
    derG_p_i = [zeros(3,3*length(ind_cont)) zeros(3,2*length(ind_slip));derG21_pi derG22_pi];
end