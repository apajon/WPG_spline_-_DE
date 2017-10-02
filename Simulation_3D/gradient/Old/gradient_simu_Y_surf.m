function Ysurf = gradient_simu_Y_surf(sole,ind_cont,Rini,Fc_mc3)
    t2 = [0 1 0];
    t1 = [1 0 0];
    Ylow = zeros(3,3*sole.nFreeSurf);
    for i=1:length(ind_cont)
        tA = t1 * Rini * Fc_mc3(3*i);
        tB = -t2 * Rini * Fc_mc3(3*i);
        tC = -t1 * Rini * Fc_mc3((3*i)-1) + t2 * Rini * Fc_mc3((3*i)-2);
        Ylow(:,(3*ind_cont(i)-2):3*ind_cont(i)) = [tA; tB; tC];
    end
    Ysurf = [zeros(3,3*sole.nFreeSurf);Ylow];    
end