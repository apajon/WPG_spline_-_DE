function Esurf = gradient_simu_E_surf(sole,ind_cont,ind_slip,Fc_mc,P_mc,Zdes)
    E11 = repmat(eye(3),1,sole.nFreeSurf);
    E21 = zeros(3,3*sole.nFreeSurf); 
    for i=1:length(ind_cont)
        E21(:,(3*ind_cont(i))-2:3*ind_cont(i)) = [0 0 -(P_mc(1,i)-Zdes(1));0 0 P_mc(2,i)-Zdes(2);-P_mc(2,i)+Zdes(2) P_mc(1,i)-Zdes(1) 0];
    end
    E22 = zeros(3,2*length(ind_slip));
    if ~isempty(ind_slip)
        Icont = [];
        Isurf = [];
        for i=1:length(ind_slip)
            Icont(i) = find(ind_cont==ind_slip(i));
            Isurf(i) = find((1:sole.nFreeSurf)==ind_slip(i));
        end
        Fc_mc_splip = Fc_mc(:,Icont);
        for i=1:length(ind_slip)
            E22(:,(2*i)-1:2*i) = [-Fc_mc_splip(3,i) 0;0 Fc_mc_splip(3,i);Fc_mc_splip(2,i) -Fc_mc_splip(1,i)];
        end
    end
    Esurf = [E11 zeros(3,2*length(ind_slip));E21 E22];
end