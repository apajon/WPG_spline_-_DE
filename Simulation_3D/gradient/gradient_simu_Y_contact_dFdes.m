function Y = gradient_simu_Y_contact_dFdes(dPc_q_1_dFdes,contAngle,ind_cont,Fc_mc3,ind)
    t1 = [1 0 0];    
    t2 = [0 1 0];    
    Y = zeros(6,3);
    if contAngle==ind
        Y(1:3,1:3) = eye(3);
    else
        for i=1:3
            nc = length(ind_cont);
            At1 = reshape(dPc_q_1_dFdes(:,i),3,nc)*Fc_mc3(1:3:3*nc);
            At2 = reshape(dPc_q_1_dFdes(:,i),3,nc)*Fc_mc3(2:3:3*nc);
            An = reshape(dPc_q_1_dFdes(:,i),3,nc)*Fc_mc3(3:3:3*nc);
            Y(:,i) = [zeros(3,1);[t1*An; -t2*An; -t1*At2+t2*At1]];
        end
    end
end