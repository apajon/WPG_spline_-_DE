function [costTorqueAnkle,cost_omega,der_costTorque_dwpg,domega_dwpg,angle_theta,angle_phi,angle_psi,dangle_theta_dwpg,dangle_phi_dwpg,dangle_psi_dwpg] = simulation(sole,param_sopt,simu,gradS,comp_grad,comp_grad_constr,gradpzmp_upd,gradfzmp_upd)
    friction = param_sopt.friction;
    Fdes = param_sopt.Fdes;
    Zdes = param_sopt.Zdes;

    sole.prep_stressVonMises(); % Preparation of VonMises' Stress computations

    % Variable Initialisation
    simu.P = [];
    simu.no_conv = [];
    P0 = sole.coor';
    simu.PfreeSurf = P0(:,sole.nodesFreeSurf);
    simu.PfreeSurf3 = reshape(simu.PfreeSurf,3*sole.nFreeSurf,1);
    % if Zdes(2,1)>0
    %     simu.displ = 1.0e-04 *[0.200672417450334;0.412701716712745;-0.813000335911768];
    %     simu.angleact = [0.002718435798210;-0.001201953074406;-0.000033898454071];
    % else
    %     simu.displ = 1.0e-04 *[0.200672417450334;0.412701716712745;-0.813000335911768];
    %     simu.angleact = [0.002718435798210;-0.001201953074406;-0.000033898454071];
    % end

    dist = zeros(size(param_sopt.MatInt,1)/2,size(param_sopt.MatInt,2));
    for i=1:(size(param_sopt.MatInt,1)/2)
        for j=1:size(param_sopt.MatInt,2)
            dist(i,j) = sqrt((Zdes(1,1)-param_sopt.MatInt(2*i-1,j))^2+(Zdes(2,1)-param_sopt.MatInt(2*i,j))^2) ;
        end
    end
    [minDist,ind] = min(dist(:));
    [m,n] = ind2sub(size(dist),ind);
    % MatInt(2*m-1:2*m,n)
    simu.displ = param_sopt.cell_displ_ini{m,n};
    simu.angleact = param_sopt.cell_angleact_ini{m,n};

    simu.displ_first = [0;0;0];
    simu.psi_first = 0;

    no_conv_tot = [];
    der_costTorque_dwpg_tot = zeros(size(gradfzmp_upd{1},2),size(Fdes,2));
    cost_tot = zeros(size(Fdes,2),1);
    K_cart_vec = cell(size(Fdes,2),1);
    cost = 0;
    cont_no_conv = 0;
    simu.no_conv = [];
    pAnkle = zeros(3,size(Fdes,2));
    torqueAnkle = zeros(3,size(Fdes,2));
    angle = zeros(3,size(Fdes,2));
    angle_theta = zeros(1,size(Fdes,2));
    angle_phi = zeros(1,size(Fdes,2));
    angle_psi = zeros(1,size(Fdes,2));
    costTorqueAnkle = 0;
    dY_dwpg_tot = cell(size(Fdes,2),1);
    dangle_theta_dwpg = zeros(size(gradfzmp_upd{1},2),size(Fdes,2));
    dangle_phi_dwpg = zeros(size(gradfzmp_upd{1},2),size(Fdes,2));
    dangle_psi_dwpg = zeros(size(gradfzmp_upd{1},2),size(Fdes,2));
    for s=1:size(Fdes,2)
        simu.Fdes = Fdes(:,s);
        simu.Zdes = Zdes(:,s);
        simu = contacts(sole,simu,friction,s);
        pAnkle(:,s) = simu.displ + simu.R * sole.pos_ankle_loc;
        angle(:,s) = simu.angleact;
        angle_theta(1,s) = angle(1,s);
        angle_phi(1,s) = angle(2,s);
        angle_psi(1,s) = angle(3,s);
        torqueAnkle(1,s) = simu.Ftot(2,1) * pAnkle(3,s) - simu.Ftot(3,1) * (pAnkle(2,s) - simu.Z(2,1));
        torqueAnkle(2,s) = simu.Ftot(3,1) * (pAnkle(1,s) - simu.Z(1,1)) - simu.Ftot(1,1) * pAnkle(3,s);
        torqueAnkle(3,s) = simu.Ftot(1,1) * (pAnkle(2,s) - simu.Z(2,1)) - simu.Ftot(2,1) * (pAnkle(1,s) - simu.Z(1,1));
        sole.stressVonMises(simu.dPabstot);
        plotsole(2,sole.elements_surf,simu.P,sole.stressVM,simu.Pc,simu.Fc,simu.Z,simu.Ftot,-37.5,30,pAnkle(:,s));
        if comp_grad
            gradS.updContVar(sole,simu,friction,s);
            gradZdesx = gradpzmp_upd{1}(s,:);
            gradZdesy = gradpzmp_upd{2}(s,:);
            gradFdesx = gradfzmp_upd{1}(s,:);
            gradFdesy = gradfzmp_upd{2}(s,:);
            gradFdesz = gradfzmp_upd{3}(s,:);
            gradS.derShape_dwpg(sole,simu,friction,gradZdesx,gradZdesy,gradFdesx,gradFdesy,gradFdesz,s);
            dpAnkle_dwpgi = zeros(3,size(gradfzmp_upd{1},2));
            for i=1:size(gradfzmp_upd{1},2)
                dpAnkle_dwpgi(:,i) = [gradS.dR_dtheta * (sole.pos_ankle_loc), gradS.dR_dphi * (sole.pos_ankle_loc), gradS.dR_dpsi * (sole.pos_ankle_loc)]*gradS.dY_dwpg(:,i) + gradS.dOl_dwpg(:,i);
                der_costTorque_dwpg_tot(i,s) = 2*((gradFdesy(1,i) * pAnkle(3,s) + simu.Ftot(2,1) * dpAnkle_dwpgi(3,i) - gradFdesz(1,i) *(pAnkle(2,s)-simu.Z(2,1)) - simu.Ftot(3,1) * dpAnkle_dwpgi(2,i) + simu.Ftot(3,1) * gradFdesy(1,i))*torqueAnkle(1) + ...
                                               (gradFdesz(1,i) * (pAnkle(1,s)-simu.Z(1,1)) + simu.Ftot(3,1) * dpAnkle_dwpgi(1,i) - simu.Ftot(3,1) * gradZdesx(1,i) - gradFdesx(1,i) * pAnkle(3,s) - simu.Ftot(1,1) * dpAnkle_dwpgi(3,i))*torqueAnkle(2)+...
                                               (gradFdesx(1,i) * (pAnkle(2,s)-simu.Z(2,1)) + simu.Ftot(1,1) * dpAnkle_dwpgi(2,i) - simu.Ftot(1,1) * gradZdesy(1,i) - gradFdesy(1,i) * (pAnkle(1,s)-simu.Z(1)) + simu.Ftot(2,1) * dpAnkle_dwpgi(1,i) - simu.Ftot(2,1) * gradZdesx(1,i))*torqueAnkle(3));
            end
            dY_dwpg_tot{s} = gradS.dY_dwpg';
            if comp_grad_constr
                dangle_theta_dwpg(:,s) = sign(simu.angleact(1)) *  gradS.dY_dwpg(1,:)';
                dangle_phi_dwpg(:,s) = sign(simu.angleact(2)) * gradS.dY_dwpg(2,:)';
                dangle_psi_dwpg(:,s) = sign(simu.angleact(3)) * gradS.dY_dwpg(3,:)';
            end            
        end

        if ~isempty(simu.no_conv)
            cont_no_conv = cont_no_conv + 1;
            no_conv_tot(cont_no_conv) = simu.no_conv;
            simu.no_conv=[];
            angle_theta(1,s) = degtorad(120);
            angle_phi(1,s) = degtorad(120);
            angle_psi(1,s) = degtorad(120);
        end
        costTorqueAnkle = costTorqueAnkle + norm(torqueAnkle)^2;
        K_cart_vec{s} = simu.Kcart;    
    end
    scaling_CTA = 3*10^6;
    costTorqueAnkle = costTorqueAnkle/scaling_CTA;

    omegax = (1/0.005) * diff(angle(1,:));
    omegay = (1/0.005) * diff(angle(2,:));
    omegaz = (1/0.005) * diff(angle(3,:));
    domega_dwpg_tot = cell(size(Fdes,2)-1,1);
    cost_omega = 0;
    domega_dwpg = zeros(size(gradfzmp_upd{1},2),1);
    for i=1:length(omegax)
        norm_omega = norm([omegax(1,i);omegay(1,i);omegaz(1,i)]);
        cost_omega = cost_omega + norm_omega^2;
        if comp_grad
            domega_dwpg_tot{i} = (1/0.005) * (dY_dwpg_tot{i+1}-dY_dwpg_tot{i});
            for j=1:size(gradfzmp_upd{1},2)
                domega_dwpg(j) = domega_dwpg(j) + 2 * (domega_dwpg_tot{i}(j,1) * omegax(1,i) + domega_dwpg_tot{i}(j,2) * omegay(1,i) + ...
                    domega_dwpg_tot{i}(j,3) * omegaz(1,i));
            end
        end
    end

    der_costTorque_dwpg = sum(der_costTorque_dwpg_tot,2);

    save 'resultsOpt/K_cart_ini.mat' K_cart_vec
end