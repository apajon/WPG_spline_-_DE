function dPs_dFdes = derPs_dFdes(sole,ind_cont3,Ws,dOl_dFdes,dY_dFdes,dFc_dFdes,b)
    W_Fc = Ws(:,ind_cont3)*dFc_dFdes;
    dPs_dFdes = repmat(dOl_dFdes,sole.nFreeSurf,1) + b * dY_dFdes + W_Fc;
end