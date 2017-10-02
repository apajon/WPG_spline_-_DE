function derG_Fdes = compute_derG_Fdes(dPc_dFdes,dFc_dFdes,ind_cont,ind_slip)
    derG21_Fdes = zeros(3,3*length(ind_cont));
    for i=1:length(ind_cont)
        derG21_Fdes(:,3*i-2:3*i) = [0 0 dPc_dFdes(3*i-1);0 0 -dPc_dFdes(3*i-2); -dPc_dFdes(3*i-1) dPc_dFdes(3*i-2) 0]; 
    end
    derG22_Fdes = [];
    if ~isempty(ind_slip)      
        Islip = zeros(length(ind_slip),1);
        for i=1:length(ind_slip)
            Islip(i) = find(ind_cont==ind_slip(i));
        end
        for i=1:length(ind_slip)
            derG22_Fdes(:,2*i-1:2*i) = [0 dFc_dFdes(3*Islip(i));dFc_dFdes(3*Islip(i)) 0; dFc_dFdes(3*Islip(i)-1) -dFc_dFdes(3*Islip(i)-2)]; 
        end    
    end
    derG_Fdes = [zeros(3,3*length(ind_cont)) zeros(3,2*length(ind_slip));derG21_Fdes derG22_Fdes];
end