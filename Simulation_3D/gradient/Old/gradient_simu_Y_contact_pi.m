function Y = gradient_simu_Y_contact_pi(ind_cont,Rini,Fc_mc3,dPFreec_dp_i,contAngle,dPcprevious_dp_i)
    t1 = [1 0 0];  
    t2 = [0 1 0];
    At1 = zeros(3,1);
    At2 = zeros(3,1);
    An = zeros(3,1);
    for j=1:length(ind_cont)
        if contAngle==1
            At1 = At1 + Rini * dPFreec_dp_i(3*j-2:3*j) * Fc_mc3(3*j-2);
            At2 = At2 + Rini * dPFreec_dp_i(3*j-2:3*j) * Fc_mc3(3*j-1);
            An = An + Rini * dPFreec_dp_i(3*j-2:3*j) * Fc_mc3(3*j);
        else
            At1 = At1 + dPcprevious_dp_i(3*j-2:3*j) * Fc_mc3(3*j-2);
            At2 = At2 + dPcprevious_dp_i(3*j-2:3*j) * Fc_mc3(3*j-1);
            An = An + dPcprevious_dp_i(3*j-2:3*j) * Fc_mc3(3*j);
        end
    end
    Y = [zeros(3,1);[t1*An; -t2*An; -t1*At2+t2*At1]];
end