function dPs_dZdes = derPs_dZdes(sole,ind_cont3,Ws,dOl_dZdes,dY_dZdes,dFc_dZdes,b)
    W_Fc = Ws(:,ind_cont3)*dFc_dZdes;
    dPs_dZdes = repmat(dOl_dZdes,sole.nFreeSurf,1) + b * dY_dZdes + W_Fc;
end