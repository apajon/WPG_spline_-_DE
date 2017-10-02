function derG_dZdes = compute_derG_Zdes(dPc_dZdes,dFc_dZdes,ind_cont,ind_slip,index_Zdes,ind,contAngle)
    derG21_dZdes = zeros(3,3*length(ind_cont));
    if ind==contAngle
        if index_Zdes==1
            for i=1:length(ind_cont)
                derG21_dZdes(:,3*i-2:3*i) = [0 0 dPc_dZdes(3*i-1);0 0 -dPc_dZdes(3*i-2)+1; -dPc_dZdes(3*i-1) dPc_dZdes(3*i-2)-1 0]; 
            end
        else
            for i=1:length(ind_cont)
                derG21_dZdes(:,3*i-2:3*i) = [0 0 dPc_dZdes(3*i-1)-1;0 0 -dPc_dZdes(3*i-2); -dPc_dZdes(3*i-1)+1 dPc_dZdes(3*i-2) 0]; 
            end        
        end
    else
        for i=1:length(ind_cont)
            derG21_dZdes(:,3*i-2:3*i) = [0 0 dPc_dZdes(3*i-1);0 0 -dPc_dZdes(3*i-2); -dPc_dZdes(3*i-1) dPc_dZdes(3*i-2) 0]; 
        end
    end
    derG22_dZdes = [];
    if ~isempty(ind_slip)      
        Islip = zeros(length(ind_slip),1);
        for i=1:length(ind_slip)
            Islip(i) = find(ind_cont==ind_slip(i));
        end
        for i=1:length(ind_slip)
            derG22_dZdes(:,2*i-1:2*i) = [0 dFc_dZdes(3*Islip(i));dFc_dZdes(3*Islip(i)) 0; dFc_dZdes(3*Islip(i)-1) -dFc_dZdes(3*Islip(i)-2)]; 
        end    
    end    
    derG_dZdes = [zeros(3,3*length(ind_cont)) zeros(3,2*length(ind_slip));derG21_dZdes derG22_dZdes];
end