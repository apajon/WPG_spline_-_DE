function Y = gradient_simu_Y_contact_dZdes(dPc_q_1_dZdes,contAngle,ind_cont,Fc_mc3,ind)
    t1 = [1 0 0];    
    t2 = [0 1 0];    
    Y = zeros(6,2);
    if contAngle==ind
        Ftotn = sum(Fc_mc3(3:3:end));
        Ftott2 = sum(Fc_mc3(2:3:end));
        Ftott1 = sum(Fc_mc3(1:3:end));
        Y(4,1)=-Ftotn;
        Y(6,1)=Ftott2;
        Y(5,2)=Ftotn;
        Y(6,2)=-Ftott1;    
    else
        for i=1:2
            nc = length(ind_cont);
            At1 = reshape(dPc_q_1_dZdes(:,i),3,nc)*Fc_mc3(1:3:3*nc);
            At2 = reshape(dPc_q_1_dZdes(:,i),3,nc)*Fc_mc3(2:3:3*nc);
            An = reshape(dPc_q_1_dZdes(:,i),3,nc)*Fc_mc3(3:3:3*nc);
            Y(:,i) = [zeros(3,1);[t1*An; -t2*An; -t1*At2+t2*At1]];
        end
    end
end