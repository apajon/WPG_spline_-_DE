function H_tilde = gradient_simu_H_tilde_contact(angleact,dY_dP,ind_cont)
    [dR_dtheta,dR_dphi,dR_dpsi] = dRot(angleact);
    dR_dtheta_block_c = [];
%     for j=1:length(ind_cont)
%         dR_dtheta_block_c = blkdiag(dR_dtheta_block_c,dR_dtheta);
%     end        
%     dR_dphi_block_c = [];
%     for j=1:length(ind_cont)
%         dR_dphi_block_c = blkdiag(dR_dphi_block_c,dR_dphi);
%     end
%     dR_dpsi_block_c = [];
%     for j=1:length(ind_cont)
%         dR_dpsi_block_c = blkdiag(dR_dpsi_block_c,dR_dpsi);
%     end
    for i=1:length(ind_cont)
        H_tilde = [dR_dtheta_block_c * dY_dP(1,:)' dR_dphi_block_c * dY_dP(2,:)' dR_dpsi_block_c * dY_dP(3,:)'];
    end
end