function [Bw,Fc_mc_hat] = computeBw(Fc_mc3,ind_cont,ind_slip,R,Wc,Pfree_c)
    delta_mc3 = Wc * Fc_mc3;
    Fc_mc_hat = zeros(3*length(ind_cont),3);
    delta_mc_hat = zeros(3*length(ind_cont),3);
 
    for i=1:length(ind_cont)
        Fc_mc_hat((3*i)-2:3*i,:) = skew(Fc_mc3(3*i-2:3*i));
        delta_mc_hat((3*i)-2:3*i,:) = skew(delta_mc3(3*i-2:3*i));
    end
    Bw12 = Wc * Fc_mc_hat - delta_mc_hat;
    for i=1:length(ind_cont)
        Bw12(3*i-2:3*i,:) = Bw12(3*i-2:3*i,:)-skew(R*Pfree_c(:,i));
    end
    Bw11 = repmat(eye(3),length(ind_cont),1);
    Bw = [Bw11 Bw12;zeros(2*length(ind_slip),3) zeros(2*length(ind_slip),3)];
end