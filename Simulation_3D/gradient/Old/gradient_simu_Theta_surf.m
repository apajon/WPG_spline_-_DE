function Theta_mat_surf = gradient_simu_Theta_surf(sole,der_C_s,ind_slip,FSurf3,RBlock_s,Rini_block_s)
    der_W_s = zeros(3*sole.nFreeSurf,3*sole.nFreeSurf);
    Rini_block_s_t = Rini_block_s;
    Rini_block_s_t(3:3:end,:)=0;
    %Rini_block_c_t(:,3:3:end)=0;
    for i=1:3*sole.nFreeSurf
        der_W_s(:,i) = RBlock_s * der_C_s{i} * RBlock_s' * FSurf3;
    end   
    Theta_mat_surf = [der_W_s+ RBlock_s - Rini_block_s_t;zeros(2*length(ind_slip),3*sole.nFreeSurf)];
end