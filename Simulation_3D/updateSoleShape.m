function [costTorqueAnkle,cost_omega,der_costTorque_dwpg,ddotY_dwpg,angle_theta,angle_phi,angle_psi,dangle_theta_dwpg,dangle_phi_dwpg,dangle_psi_dwpg] = updateSoleShape(sole,param_sopt,comp_grad,comp_grad_constr,simu,gradS,gradpzmp_upd,gradfzmp_upd)

% Start Simulation
[costTorqueAnkle,cost_omega,der_costTorque_dwpg,ddotY_dwpg,angle_theta,angle_phi,angle_psi,dangle_theta_dwpg,dangle_phi_dwpg,dangle_psi_dwpg] = simulation(sole,param_sopt,simu,gradS,comp_grad,comp_grad_constr,gradpzmp_upd,gradfzmp_upd);
end