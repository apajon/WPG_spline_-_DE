function cost = testder(p_ini_v,sole,Fdes,Zdes,friction,spl,move_dirichlet,w)
% Create structure for the shape optimization parameters
simu = Simu(); %create Simulation object 

[sole.coor,~] = deformation_moveDiri(sole,p_ini_v,spl,move_dirichlet);

sole.stiffness();

sole.prep_stressVonMises(); % Preparation of VonMises' Stress computations

% Variable Initialisation
simu.P = [];
simu.no_conv = [];
P0 = sole.coor';
simu.PfreeSurf = P0(:,sole.nodesFreeSurf);
simu.PfreeSurf3 = reshape(simu.PfreeSurf,3*sole.nFreeSurf,1);
simu.displ = [0.00000262899049116;-0.00000406461222954;-0.0001155602271939];
simu.angleact = [-0.126916077961588;-0.118047287866647;-0.000913850459974];
simu.displ_first = [0;0;0];
simu.psi_first = 0;
cost_lin = 0;
cost_rot = 0;

for s=1:2
    simu.Fdes = Fdes(:,s);
    simu.Zdes = Zdes(:,s);
    [simu,~] = contacts(sole,simu,friction,s);
    
    sole.stressVonMises(simu.dPabstot);
    plotsole(2,sole.elements_surf,simu.P,sole.stressVM,simu.Pc,simu.Fc,simu.Z,simu.Ftot,-37.5,30);


    cost_rot = cost_rot + sqrt(min(eig(simu.Kcart(4:6,4:6)'*simu.Kcart(4:6,4:6)))); % maximize the minimum stiffness. Be carefull, this is not differentiable when changing of minimum eigenvalue !!!!!!
%     if comp_grad
%         if s==1
%            [der_Kctr_dp,der_Kcrot_dp,dPs_dp] = compute_gradient_step(sole,size(dPFree_dp,2),Kcart,displ,angleact,Zdes,friction,PfreeSurf,PSurf3,PabsOld,Fc_mc3,sole.Cs,dPFree_dp,ind_cont,ind_slip,displ_first,psi_first,s,w,B_alg,[]);
%         else
%             [der_Kctr_dp,der_Kcrot_dp,dPs_dp] = compute_gradient_step(sole,size(dPFree_dp,2),Kcart,displ,angleact,Zdes,friction,PfreeSurf,PSurf3,PabsOld,Fc_mc3,sole.Cs,dPFree_dp,ind_cont,ind_slip,displ_first,psi_first,s,w,B_alg,dPs_dp);
%         end
%     end
% 
    if s<(size(Fdes,2)/10) %minimize linear stiffness for the first 10% of the step
        cost_lin = cost_lin + sqrt(max(eig(simu.Kcart(1:3,1:3)'*simu.Kcart(1:3,1:3))));
%         if comp_grad
%             der_Kctr_dptot = der_Kctr_dptot + der_Kctr_dp;
%         end
    end
%     
%     if comp_grad
%         der_Kcrot_dptot = der_Kcrot_dptot + der_Kcrot_dp;
%     end    
end
cost = cost_lin - w * cost_rot;
%Ps3 = simu.PSurf3;
end