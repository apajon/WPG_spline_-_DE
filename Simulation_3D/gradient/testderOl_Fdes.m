function cost = testderOl_Fdes(Fdesx,Fdesy,Fdesz,sole,Fdes,Zdes,Zdesx,Zdesy,fric,angleact,displ,SimAna,w)
%[sole.coor,dPFree_dp] = deformation_moveDiri(sole,p_ini_v,spl,move_dirichlet);
%sole.coor = coorini;
sole.stiffness();
sole.stiffnessSurface();
%w=1;
D=3;
cost = 0;
for i=1:10
    contAngle = i;
%     a=load('Pg.mat');
%     Pg = a.Pg;
%     b = load('dPs_dp.mat');
%     dPs_dp = b.dPs_dp;
%     c = load('displ.mat');
%     displ = c.displ;
%     d = load('angleact.mat');
%     angleact = d.angleact;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    m = sole.nFreeSurf;
    Cs = sole.Cs;
    %no_conv = [];
    P0 = sole.coor';
    if contAngle==1
        Pfree = P0(:,sole.nodesFreeSurf);
        PabsOld(1,:) = P0(1,sole.nodesFreeSurf);
        PabsOld(2,:) = P0(2,sole.nodesFreeSurf);
        PabsOld(3,:) = P0(3,sole.nodesFreeSurf);    
        FSurfNew3 = zeros(D*m,1);
        displ_first= [0;0;0];
        psi_first = 0;
    else
        Pfree = P0(:,sole.nodesFreeSurf);
        PabsOld = Pg(sole.nodesFreeSurf,:)';
        %FcNew = FcFreeSurf'; 
        %FSurfNew3 = reshape(FcNew,D*m,1);
        FSurfNew3 = zeros(D*m,1);
        displ_first= [0;0;0];
        psi_first = 0;
    end
    if i==SimAna
        FtotZMPdes = [Fdesx;Fdesy;Fdesz;Zdesx;Zdesy;0]; 
    else
        FtotZMPdes = [Fdes(:,contAngle);Zdes(:,contAngle);0]; 
    end
    [displ,angleact,FtotZMP,Fc_mc3_out,PSurf3,ind_cont_out,ind_slip_out,Kcart,Bomega_out,J2,A_out,B_out] = GaussFtotZMP(contAngle,fric,m,displ,angleact,FtotZMPdes,FSurfNew3,Pfree,PabsOld,Cs,displ_first);
    theta = angleact(1);
    phi = angleact(2);
    psi = angleact(3);
    %R = Rot(theta,phi,psi);
    if ind_cont_out(1)==0
        a = find(ind_cont_out==0,2);
        Fc_mc3 = Fc_mc3_out(1:D*(a(2)-1));
        ind_cont = (sort(ind_cont_out(1:(a(2)-1))+1))';      
    else
        a = find(ind_cont_out==0,1);
        Fc_mc3 = Fc_mc3_out(1:D*(a-1));
        ind_cont = (sort(ind_cont_out(1:(a-1))+1))';              
    end
    if norm(ind_slip_out)>0
        if ind_slip_out(1)==0
            b = find(ind_slip_out==0,2);
            ind_slip = (sort(ind_slip_out(1:(b(2)-1))+1))';  
        else
            b = find(ind_slip_out==0,1);
            ind_slip = (sort(ind_slip_out(1:(b-1))+1))';      
        end
    else
        ind_slip = [];
    end
    ind_cont3 = sort(([D*ind_cont-2; D*ind_cont-1; D*ind_cont]));
    ind_cont3 = reshape(ind_cont3,D*length(ind_cont),1); 
    Pfree_c = Pfree(:,ind_cont);
    Ftot = FtotZMP(1:3);
    Z = FtotZMP(4:5);
    theta = angleact(1);
    phi = angleact(2);
    psi = angleact(3);
    R = Rot(theta,phi,psi);
    Rini = Rot(theta,phi,psi_first);
    cost_rot = sqrt(min(eig(Kcart(4:6,4:6)'*Kcart(4:6,4:6))));
    cost_lin = sqrt(max(eig(Kcart(1:3,1:3)'*Kcart(1:3,1:3))));
    cost = cost + (cost_lin - w * cost_rot);
    Fc_mc = reshape(Fc_mc3,3,length(ind_cont));
    P_mc3 = PSurf3(ind_cont3);
    P_mc = reshape(P_mc3,3,length(ind_cont));   
    nf = sole.nFreeSurf;
    tmp = reshape(R*reshape(Cs,3,3*nf*nf),3*nf,3*nf);
    tmp2 = tmp';
    Ws = reshape(R*reshape(tmp2,3,3*nf*nf),3*nf,3*nf);
    Cc = Cs(ind_cont3,ind_cont3);
    Wc = Ws(ind_cont3,ind_cont3);
    PabsOld_mc = PabsOld(:,ind_cont);
%     J2_alg=J2(1:(6*6));
%     J2_alg=reshape(J2_alg,6,6);
%     A_alg = A_out(1:(D*length(ind_cont)+D*length(ind_slip))*(D*length(ind_cont)+D*length(ind_slip)));
%     A_alg = reshape(A_alg,D*length(ind_cont)+D*length(ind_slip),D*length(ind_cont)+D*length(ind_slip));
%     Pfree_s = Pfree;
%     if contAngle==1
%         [der_Kctr_dp,der_Kcrot_dp,dPs_dp] = compute_gradient_step(sole,size(dPFree_dp,2),displ,angleact,Zdes,fric,Pfree,PSurf3,PabsOld,Fc_mc3,Cs,dPFree_dp,ind_cont,ind_slip,displ_first,psi_first,contAngle,w,[]);
%     else
%         [der_Kctr_dp,der_Kcrot_dp,dPs_dp] = compute_gradient_step(sole,size(dPFree_dp,2),displ,angleact,Zdes,fric,Pfree,PSurf3,PabsOld,Fc_mc3,Cs,dPFree_dp,ind_cont,ind_slip,displ_first,psi_first,contAngle,w,dPs_dp);
%     end
%     der_cost_dp = der_Kctr_dp - der_Kcrot_dp;
%     %angleact = testderY_Pfree_i(sole,coorini);
%     %f = @(x) testderY_Pfree_i(x,sole,Zdes(:,1),Fdes(:,1),fric,angleact,displ);
%     f = @(x) testderOl_Pfree_i(x,sole,Zdes(:,contAngle),Fdes(:,contAngle),fric,angleact,displ,spl,move_dirichlet,contAngle);
%     J5 = diff5points(f,1,1e-4,p_ini_v);
%     (dPc_dPfree(:,4)-J{2})./dPc_dPfree(:,4)
    %H = 10.^[-10:0.1:-2];
%     figure;
%     for i=1:141
%         %y = []; 
% %         for h=H
%             %J5 = diff5points(f,i,1e-4,p_ini_v); 
%             J = der_cost_dp(i,1); 
%             y = (J-J5{i})/J;
% %         end; 
%          
%         stem(i,y);
%         hold on
%     end;
%     H = 10.^[-10:0.1:-2];
%     figure;
%     for i=121
%         %y = []; 
%         for h=H
%             J5 = diff5points(f,[121,122],h,p_ini_v);
%             J = der_cost_dp(i,1); 
%             y = (J-J5{i})/J;
%             stem(h,y);
%             hold on            
%         end; 
%          
% 
%   end
    Xi = computeXi(displ,Zdes);
    A = gradient_simu_A_contact(P_mc,Fc_mc,ind_cont,Wc,ind_slip,PabsOld_mc,fric,contAngle,Rini,displ_first);
    G = computeG(ind_cont,ind_slip,Fc_mc,P_mc,Zdes);
    [Bw,Fc_mc_hat] = computeBw(Fc_mc3,ind_cont,ind_slip,R,Wc,Pfree_c);
    %% Find the displacement of internal nodes and position of all nodes
    FcSurf3 = zeros(D*sole.nFreeSurf,1);
    FcSurf3(ind_cont3) = Fc_mc3;
    FcFreeSurf = reshape(FcSurf3,D,m)';
    Fcloc = R'*FcFreeSurf';
    Fcloc3 = reshape(Fcloc,D*m,1);
    dPlocSurf3 = sole.Cs * Fcloc3;
    %Displacement of all nodes
    dPInt3 = -(sole.m_invKii_Kis) * dPlocSurf3;
    dPloc3 = zeros(D*sole.nTot,1);
    dPloc3(sole.nodesFreeSurf3) = dPlocSurf3;
    dPloc3(sole.nodesInt3) = dPInt3;
    dPloc = reshape(dPloc3,D,sole.nTot)';
    dPabstot = R*dPloc';
    displacement = repmat(displ', sole.nTot, 1);
    %PosRot = R*P0 + center;
    PosRot = R*sole.coor';
    Pg = PosRot' + dPabstot' + displacement;
    %Pg(sole.nodesFreeSurf,:)=PSurf';
    Fc = FcFreeSurf(ind_cont,:);
    Pc = Pg(sole.nodesFreeSurf(ind_cont),:);
    
    stressVM0 = zeros(sole.nTot,1); % VonMises' Stress    
    plotsole(2,sole.elements_surf,Pg,stressVM0,Pc,Fc,Z,Ftot,-127.5,30);

end
% theta = angleact(1);
% phi = angleact(2);
% psi = angleact(3);
% R = Rot(theta,phi,psi);
% Rini = Rot(theta,phi,0); 
% P_mc = reshape(P_mc3,3,length(ind_cont));
% Fc_mc = reshape(Fc_mc3,3,length(ind_cont));
% Wc = zeros(3*length(ind_cont),3*length(ind_cont));
% for j = 1:length(ind_cont)
%     for k = 1:length(ind_cont)
%         Wc((3*j)-2:3*j,(3*k)-2:3*k) = R * Cs((3*ind_cont(j))-2:3*ind_cont(j),(3*ind_cont(k))-2:3*ind_cont(k)) * R';
%     end
% end
% Pfree_c = Pfree(:,ind_cont);
% A = gradient_simu_A_contact(P_mc,Fc_mc,ind_cont,Wc,ind_slip,Pfree_c,fric,contAngle,Rini,displ_first);
% Bw = computeBw(Fc_mc3,ind_cont,ind_slip,R,Wc,Pfree_c); 
% G = computeG(ind_cont,ind_slip,Fc_mc,P_mc,Zdes);
end