function Y = gradient_simu_Y_contact(ind_cont,Rini,Fc_mc3,dPFreec_dp,lp,contAngle,dPcq_1_dp)
    t1 = [1 0 0];    
    t2 = [0 1 0];   
    Y = zeros(6,lp);
    if contAngle==1
        for i=1:lp
            At1 = zeros(3,1);
            At2 = zeros(3,1);
            An = zeros(3,1);
            for j=1:length(ind_cont)
                At1 = At1 + Rini * dPFreec_dp(3*j-2:3*j,i) * Fc_mc3(3*j-2);
                At2 = At2 + Rini * dPFreec_dp(3*j-2:3*j,i) * Fc_mc3(3*j-1);
                An = An + Rini * dPFreec_dp(3*j-2:3*j,i) * Fc_mc3(3*j);
            end
            Y(:,i) = [zeros(3,1);[t1*An; -t2*An; -t1*At2+t2*At1]];
        end
    else 
        for i=1:lp
            %for j=1:length(ind_cont)
            %    At1 = At1 + dPcq_1_dp(3*j-2:3*j,i) * Fc_mc3(3*j-2);
            %    At2 = At2 + dPcq_1_dp(3*j-2:3*j,i) * Fc_mc3(3*j-1);
            %    An = An + dPcq_1_dp(3*j-2:3*j,i) * Fc_mc3(3*j);
            %end
            nc = length(ind_cont);
            At1 = reshape(dPcq_1_dp(:,i),3,nc)*Fc_mc3(1:3:3*nc);
            At2 = reshape(dPcq_1_dp(:,i),3,nc)*Fc_mc3(2:3:3*nc);
            An = reshape(dPcq_1_dp(:,i),3,nc)*Fc_mc3(3:3:3*nc);
            Y(:,i) = [zeros(3,1);[t1*An; -t2*An; -t1*At2+t2*At1]];
        end
    end
end