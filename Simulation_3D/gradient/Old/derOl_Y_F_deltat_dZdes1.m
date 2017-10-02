function  [dOl_dZdes,dY_dZdes,dFc_dZdes,ddeltat_dZdes] = derOl_Y_F_deltat_dZdes1(ind_cont,EA_1B,A_1B,Fc_mc3)
    Ftotn = sum(Fc_mc3(3:3:end));
    Ftott2 = sum(Fc_mc3(2:3:end));
    Ftott1 = sum(Fc_mc3(1:3:end));
    Y = zeros(6,2);
    Y(4,1)=-Ftotn;
    Y(6,1)=Ftott2;
    Y(5,2)=Ftotn;
    Y(6,2)=-Ftott1;    
    dOl_dY = EA_1B \ Y;
    dF_ddelta_dZdes = A_1B * dOl_dY;
    dFc_dZdes = dF_ddelta_dZdes(1:3*length(ind_cont),:);
    ddeltat_dZdes = dF_ddelta_dZdes(3*length(ind_cont)+1:end,:);

    dOl_dZdes = dOl_dY(1:3,:);
    dY_dZdes = dOl_dY(4:6,:);
end