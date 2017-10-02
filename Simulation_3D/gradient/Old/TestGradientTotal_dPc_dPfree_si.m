format long
clc
clear all
close all

addpath ./input
addpath ./FEM
addpath ./iso2mesh
addpath ./geom3d/geom3d
addpath ./plane_line_intersect
addpath ./gradient

clear sole soleini
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Sole FEM                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fname = 'semelle1.msh';
pname = 'input/semelle1 L=0.23, l=0.13, e=0.03 m new centre less detailed/';
%%% foot size %%%
l=0.06;
L=0.13;
e=0.06;
Point = zeros(8,3);
Point(1,:) = [0, -l/2, 0];
Point(2,:) = [0, l/2, 0];
Point(3,:) = [0, l/2, e];
Point(4,:) = [0, -l/2, e];
Point(5,:) = [L, -l/2, 0];
Point(6,:) = [L, l/2, 0];
Point(7,:) = [L, l/2, e];
Point(8,:) = [L, -l/2, e];
sole = soleFEM_newStiff(pname,fname,l,L,e);
soleini = soleFEM_newStiff(pname,fname,l,L,e);
%%% In the gmesh format the Point are in the first positions of sole.coor
for i=1:sole.nTot
    if sole.coor(i,:)==Point(1,:);
        I_point1 = i;
    end
    if sole.coor(i,:)==Point(2,:);
        I_point2 = i;
    end
    if sole.coor(i,:)==Point(3,:);
        I_point3 = i;
    end
    if sole.coor(i,:)==Point(4,:);
        I_point4 = i;
    end
    if sole.coor(i,:)==Point(5,:);
        I_point5 = i;
    end
    if sole.coor(i,:)==Point(6,:);
        I_point6 = i;
    end
    if sole.coor(i,:)==Point(7,:);
        I_point7 = i;
    end  
    if sole.coor(i,:)==Point(8,:);
        I_point8 = i;
    end      
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           B-spline                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
spline_res = -pi/2:pi/10:pi/2;
l_spli = length(spline_res);
spl = SplineClass(sole,spline_res);
splini = SplineClass(soleini,spline_res);
%%%%%%%%%%%%%%%%%%%%%%%%%%0%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        Shape Optimization                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
type = 1;
[move_dirichlet,first_opt,trajectory,material,w,s1,s2,s3,p_ini_v] = setConfiguration(type,l_spli);

% fig1 = figure(2); clf;
% set(gcf,'PaperUnits','inches','PaperPosition',[0 0 40 30])
% mTextBox = uicontrol('style','text');
% set(mTextBox,'String','Material=Neoprene','BackgroundColor',[1 1 1],'Position',[55 2315 2000 200],'fontsize',100);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          Material                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%(%%%%%%%%%%%%%%%%%
% Elastomer
% Butyl Rubber 0.001-0.002 GPa
% Sylicon Elastomers 0.005-0.02 GPa
% Neoprene (CR) 0.0007-0.002 GPa
% Neoprene (CR)
if material == 1
    Young = 1000000;
else
    Young = 1400000;    
end

Poisson = 0.3;
sole.setMaterial(Young,Poisson);
soleini.setMaterial(Young,Poisson);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               ZMP                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Linux
% [t, ZMPx, ZMPy, COMx, COMy, Fx, Fy, Fz] = textread('/media/63DFD35C3FD8A3D0/Giappone/Code/sim3d/Simulation 3D - Desired ZMP - Position and Force - 3 Angles - meters - Linux/trajectory/exemple_trajectoire.txt', '%f %f %f %f %f %f %f %f', 'headerlines', 1);
% Windows
if exist('trajectory/exemple_trajectoire.txt','file')~=0
    % Linux
    [t, ZMPx, ZMPy, COMx, COMy, Fx, Fy, Fz] = textread('trajectory/exemple_trajectoire.txt', '%f %f %f %f %f %f %f %f', 'headerlines', 1);
elseif exist('trajectory\exemple_trajectoire.txt','file')~=0
    % Windows
    [t, ZMPx, ZMPy, COMx, COMy, Fx, Fy, Fz] = textread('trajectory\exemple_trajectoire.txt', '%f %f %f %f %f %f %f %f', 'headerlines', 1);
end

backtoankle=0.098; %from back to ankle of foot
fronttoankle=0.128; %from  front to ankle of foot
exttoankle=0.076; %from exterior to ankle of foot
inttoankle=0.054; %from interior to ankle of foot
xpankle=1.257728354176543; %x coordinate of ankle position
ypankle=-0.045000000000611; %y coordinate of ankle position
rightorleft = 1;
[Zdes,Fdes] = changeRef(ZMPx, ZMPy, Fx, Fy, Fz,backtoankle,fronttoankle,exttoankle,inttoankle,xpankle,ypankle,rightorleft,trajectory);

% print(fig1,'results/CMAME/Video/IniSlideNeTr2','-dpng','-r100')
Zdes(:,end-10:end)=[];
Zdes(:,1:10)=[];
Zdes = Zdes(:,1:1:length(Zdes));

Fdes(:,end-10:end)=[];
Fdes(:,1:10)=[];
Fdes = Fdes(:,1:1:length(Fdes));

%%% Create shape parameters struct
fric = 0.8;

%sole.coor = deformation_moveDiri(sole,p_ini_v,spl,param_sopt.move_dirichlet); % Moving Dirichlet
st1 = 1e-09;
% p_ini_v(9)=p_ini_v(9) + st1;
[sole.coor,dPtot3_dp] = deformation_moveDiri(sole,p_ini_v,spl,move_dirichlet);
dPFreeSurf3_dp = dPtot3_dp(sole.nodesFreeSurf3,:);
stressVM0 = zeros(sole.nTot,1);
plotsole(1,sole.elements_surf,sole.coor,stressVM0,[],[],[],[],-37.5,30);

% sole.coor(1,1) = sole.coor(1,1) + st1;
sole.stiffness();
sole.stiffnessSurface();
sole.derStiff();
% a = load('Pg.mat');
% Pg = a.Pg();
% b = load('FcFreeSurf.mat');
% FcFreeSurf = b.FcFreeSurf();
% a = load('displ.mat');
% displ = a.displ();
% b = load('angleact.mat');
% angleact = b.angleact();
% displ = 1.0e-03 * [0.115457331163844;0.006202241467888;-0.512942670322166];
% angleact = [0.000058485776918;-0.007589944959809;-0.000058477583626];
ABe = prep_stressVonMises(sole);
%%% Test Gradient algo with given position and orientation
a = load('displ_tot.mat');
displ_tot = a.displ_tot(); 
b = load('angleact_tot.mat');
angleact_tot = b.angleact_tot();
Pg = [];
FcFreeSurf = [];
D = 3;
for i=1:size(angleact_tot,2)
    angleact = angleact_tot(:,i);
    displ = displ_tot(:,i);
    contAngle = i;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PsoleAngle.m
    m = sole.nFreeSurf;
    Cs = sole.Cs;
    no_conv = [];
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
        FcNew = FcFreeSurf'; 
        FSurfNew3 = reshape(FcNew,D*m,1);
        displ_first= [0;0;0];
        psi_first = 0;
    end
    FcOld3 = FSurfNew3;
    FcOld = reshape(FcOld3,3,sole.nFreeSurf);
    diplini = displ;
    angleactini = angleact;
    angleactOld = angleact;
    Fdes(:,i) = Fdes(:,i)*100;
    FtotZMPdes = [Fdes(:,i);Zdes(:,i);0]; 
    [displ,angleact,FtotZMP,Fc_mc3_out,PSurf3,ind_cont_out,ind_slip_out,Kcart,Bomega_out,J2,A_out,B_out] = GaussFtotZMP(contAngle,fric,m,displ,angleact,FtotZMPdes,FSurfNew3,Pfree,PabsOld,Cs,displ_first);

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

    Ftot = FtotZMP(1:3);
    Z = FtotZMP(4:5);
    ind_cont3 = sort(([D*ind_cont-2; D*ind_cont-1; D*ind_cont]));
    ind_cont3 = reshape(ind_cont3,D*length(ind_cont),1); 
    Pfree_c = Pfree(:,ind_cont);
    Pfree_c3 = reshape(Pfree_c,3*length(ind_cont),1);
    Fc_mc = reshape(Fc_mc3,3,length(ind_cont));
    PSurf = reshape(PSurf3',3,sole.nFreeSurf);
    P_mc3 = PSurf3(ind_cont3);
    P_mc = reshape(P_mc3,3,length(ind_cont));
%      dPsurf_dPfrees(:,sole.nodesFreeSurf3(10))
    PabsOld_New = PabsOld(:,ind_cont);
    PabsOld_New3 = reshape(PabsOld_New,3*length(ind_cont),1);
    theta = angleact(1);
    phi = angleact(2);
    psi = angleact(3);
    R = Rot(theta,phi,psi);
    RBlock_c = [];
    for j=1:length(ind_cont)
        RBlock_c = blkdiag(RBlock_c,R);
    end
    RBlock_s = [];
    for j=1:sole.nFreeSurf
        RBlock_s = blkdiag(RBlock_s,R);
    end    
    psi_first = 0;
    Rini = Rot(theta,phi,psi_first);
    Rini_block_c = [];
    for j=1:length(ind_cont)
        Rini_block_c = blkdiag(Rini_block_c,Rini);
    end
    Rini_block_s = [];
    for j=1:sole.nFreeSurf
        Rini_block_s = blkdiag(Rini_block_s,Rini);
    end        
    Wc = zeros(3*length(ind_cont),3*length(ind_cont));
    for j = 1:length(ind_cont)
        for k = 1:length(ind_cont)
            Wc((3*j)-2:3*j,(3*k)-2:3*k) = R * Cs((3*ind_cont(j))-2:3*ind_cont(j),(3*ind_cont(k))-2:3*ind_cont(k)) * R';
        end
    end
    PabsOld_mc = PabsOld(:,ind_cont);
    A = gradient_simu_A_contact(P_mc,Fc_mc,ind_cont,Wc,ind_slip,PabsOld_mc,fric,contAngle,displ_first,psi_first,angleact);
    Ws = zeros(3*sole.nFreeSurf,3*sole.nFreeSurf);
    for j = 1:sole.nFreeSurf
        for k = 1:sole.nFreeSurf
            Ws((3*j)-2:3*j,(3*k)-2:3*k) = R * Cs((3*j)-2:3*j,(3*k)-2:3*k) * R';
        end
    end    
%     A_surf = gradient_simu_A_surf(sole,P_mc,Fc_mc,ind_cont,ind_cont3,Wc,ind_slip,PabsOld_mc,fric,contAngle,displ_first,psi_first,angleact);

    Ccc = Cs(ind_cont3,ind_cont3);
    B = gradient_simu_B_contact(angleact,Pfree_c3,Fc_mc3,ind_cont,Ccc,ind_slip,psi_first,contAngle);
    Pfree_s = Pfree;
    Pfree_s3 = reshape(Pfree_s,3*sole.nFreeSurf,1);
    FSurf3 = zeros(D*sole.nFreeSurf,1);
    FSurf3(ind_cont3) = Fc_mc3;
%     FSurf = reshape(FcSurf3,D,m)';
%     B_surf = gradient_simu_B_surf(sole,angleact,Pfree_s3,FSurf3,Cs,ind_slip,psi_first,contAngle);
%     A_alg = A_out(1:(D*length(ind_cont)+D*length(ind_slip))*(D*length(ind_cont)+D*length(ind_slip)));
%     A_alg = reshape(A_alg,D*length(ind_cont)+D*length(ind_slip),D*length(ind_cont)+D*length(ind_slip));
%     B_alg = B_out(1:6*(D*length(ind_cont)+D*length(ind_slip)));
%     B_alg = reshape(B_alg,D*length(ind_cont)+D*length(ind_slip),6);
%     G_alg = G_out(1:6*(D*length(ind_cont)+D*length(ind_slip)));
%     G_alg = reshape(G_alg,6,D*length(ind_cont)+D*length(ind_slip));
%     J2_alg = J2(1:(6*6));
%     J2_alg = reshape(J2_alg,6,6);
    der_C_s = cell(3*sole.nFreeSurf,1);
    KNodesFreeSurf = reshape(sole.dof(sole.nodesFreeSurf,:)',[1,3*sole.nFreeSurf]);
    for j = 1:(3*sole.nFreeSurf)      
        der_C_s{j} = sole.der_C{KNodesFreeSurf(j)}(KNodesFreeSurf,KNodesFreeSurf);
    end

    der_C_cc = cell(3*length(ind_cont),1);
    for j = 1:(3*length(ind_cont))
        der_C_cc{j} = der_C_s{ind_cont3(j)}(ind_cont3,ind_cont3);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  Ol_dPfree and dY_dPfree   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    Theta_mat = gradient_simu_Theta_contact(ind_cont,der_C_cc,ind_slip,Fc_mc3,RBlock_c,Rini_block_c);
%     Theta_mat_surf = gradient_simu_Theta_surf(sole,der_C_s,ind_slip,FSurf3,RBlock_s,Rini_block_s);
    E = gradient_simu_E_contact(ind_cont,ind_slip,Fc_mc,P_mc,Zdes);
%     E_surf = gradient_simu_E_surf(sole,ind_cont,ind_slip,Fc_mc,P_mc,Zdes);
    N_ini = gradient_simu_N_ini_contact(ind_cont,Rini,Pfree_c3,Fc_mc3,angleact,psi_first,contAngle);
%     N_ini_surf = gradient_simu_N_ini_surf(ind_cont,Rini,Pfree_c3,Fc_mc3,angleact,psi_first,contAngle);
    Y = gradient_simu_Y_contact(ind_cont,Rini,Fc_mc3);
%     Y_surf = gradient_simu_Y_surf(sole,ind_cont,Rini,Fc_mc3);
%     Nini = gradient_simu_Nini_contact(ind_cont,Rini,P_mc3,Fc_mc3,angleact,psi_first,contAngle);
    %G1 = inv(A) * Theta_mat;
%     Grad1 = inv((E / A * B) + N_ini)*((-E / A *Theta_mat) + Y);
%     %Grad1 = inv((E / A *B) + Nini)*((-E / A *Theta_mat));
%     Gc = Grad1(1:3,1);
%   -0.013923711816379
%    0.018018773252787
%   -0.072736573711783

    inv_A = inv(A);
    A_1B = inv_A * B;
    inv_A_theta = inv_A * Theta_mat;
    dOl_dY = inv((E * A_1B) + N_ini)*((-E * inv_A_theta) + Y);
    
% %     %%% surf %%%
%     inv_A_surf = inv(A_surf);
%     A_1B_surf = inv_A_surf * B_surf;
%     inv_A_theta_surf = inv_A_surf * Theta_mat_surf;
%     dOl_dY_surf = inv((E_surf * A_1B_surf) + N_ini_surf)*((-E_surf * inv_A_theta_surf) + Y_surf);    
    
    dOl_dPfree = dOl_dY(1:3,:);
    dOl_dPfree_surf = zeros(3,3*sole.nFreeSurf);
    dOl_dPfree_surf(:,ind_cont3) = dOl_dPfree;
    dY_dPfree = dOl_dY(4:6,:);
    dY_dPfree_surf = zeros(3,3*sole.nFreeSurf);
    dY_dPfree_surf(:,ind_cont3) = dY_dPfree;    
    %Gc = dOl_dY(1:3,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %    dF_dPfree and ddelta_dPfree     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    dF_ddelta = A_1B * dOl_dY + inv_A_theta;
    dF_dPfree = dF_ddelta(1:3*length(ind_cont),:);
%     dF_dPfree_surf = A_1B_surf * dOl_dY_surf + inv_A_theta_surf;
    dF_dPfree_surf = zeros(3*sole.nFreeSurf,3*sole.nFreeSurf);
    for j=1:3*length(ind_cont)
        dF_dPfree_surf(ind_cont3,ind_cont3(j)) = dF_dPfree(:,j);
    end
%     KInternalNodes = reshape(obj.dof(obj.nodesInt,:)',[1,3*obj.nInt]);
%     KNodesFreeSurf = reshape(obj.dof(obj.nodesFreeSurf,:)',[1,3*obj.nFreeSurf]);
%     KNodesFree = reshape(obj.dof(obj.nodesFree,:)',[1,3*obj.nFree]);       
    dF_dPfree_surf_int = zeros(3*sole.nFree,3*sole.nFree);
    for j=1:3*sole.nFreeSurf
        dF_dPfree_surf_int(sole.nodesFreeSurf3,sole.nodesFreeSurf3(j)) = dF_dPfree_surf(:,j);
    end    
    ddelta_dPfree = dF_ddelta(3*length(ind_cont)+1:end,:);
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%     %    dPc_dPfreec             %
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     dPc_dPfreec = zeros(3*length(ind_cont),3*length(ind_cont));
%     [dR_dtheta,dR_dphi,dR_dpsi] = dRot(angleact);
%     for g=1:length(ind_cont)
%         for j=1:length(ind_cont)
%             btheta = dR_dtheta*Pfree_c3((3*j-2):3*j);
%             bphi = dR_dphi*Pfree_c3((3*j-2):3*j);
%             bpsi = dR_dpsi*Pfree_c3((3*j-2):3*j);
%             W_Fc = zeros(3,3);
%             RdCRn = zeros(3,1);
%             RdCRt1 = zeros(3,1);
%             RdCRt2 = zeros(3,1);
%             for k=1:length(ind_cont)
%                 btheta = btheta + (dR_dtheta * Ccc((3*j-2):3*j,(3*k-2):3*k) * R' + R * Ccc((3*j-2):3*j,(3*k-2):3*k) * dR_dtheta') * Fc_mc3((3*k-2):3*k);
%                 bphi = bphi + (dR_dphi * Ccc((3*j-2):3*j,(3*k-2):3*k) * R' + R * Ccc((3*j-2):3*j,(3*k-2):3*k) * dR_dphi') * Fc_mc3((3*k-2):3*k);
%                 bpsi = bpsi + (dR_dpsi * Ccc((3*j-2):3*j,(3*k-2):3*k) * R' + R * Ccc((3*j-2):3*j,(3*k-2):3*k) * dR_dpsi') * Fc_mc3((3*k-2):3*k);
%                 W_Fc = W_Fc + Wc((3*j-2):3*j,(3*k-2):3*k) * dF_dPfree((3*k-2):3*k,(3*g-2):3*g); 
%                 RdCRt1 = RdCRt1 + R * der_C_cc{(3*g)-2}((3*j-2):3*j,(3*k-2):3*k) * R' * Fc_mc3((3*k-2):3*k);
%                 RdCRt2 = RdCRt2 + R * der_C_cc{(3*g)-1}((3*j-2):3*j,(3*k-2):3*k) * R' * Fc_mc3((3*k-2):3*k);
%                 RdCRn = RdCRn + R * der_C_cc{3*g}((3*j-2):3*j,(3*k-2):3*k) * R' * Fc_mc3((3*k-2):3*k);
%             end
%             A=zeros(3,3);
%             if j==g
%                A=R;
%             end
%             dPc_dPfreec((3*j-2):3*j,(3*g-2):3*g) = dOl_dPfree(:,(3*g-2):3*g) + [btheta bphi bpsi] * dY_dPfree(:,(3*g-2):3*g) + A + W_Fc + [RdCRt1 RdCRt2 RdCRn];
%         end
%     end
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%        dPc_dPfrees              %
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dPc_dPfrees = zeros(3*length(ind_cont),3*sole.nFreeSurf);
    [dR_dtheta,dR_dphi,dR_dpsi] = dRot(angleact);
    contj = 1;
    flag=0;
    for g=1:sole.nFreeSurf
        for j=1:length(ind_cont)
            btheta = dR_dtheta*Pfree_c3((3*j-2):3*j);
            bphi = dR_dphi*Pfree_c3((3*j-2):3*j);
            bpsi = dR_dpsi*Pfree_c3((3*j-2):3*j);
            W_Fs = zeros(3,3);
            RdCRn = zeros(3,1);
            RdCRt1 = zeros(3,1);
            RdCRt2 = zeros(3,1);
            for k=1:sole.nFreeSurf
                btheta = btheta + (dR_dtheta * Cs((3*ind_cont(j)-2):3*ind_cont(j),(3*k-2):3*k) * R' + R * Cs((3*ind_cont(j)-2):3*ind_cont(j),(3*k-2):3*k) * dR_dtheta') * FSurf3((3*k-2):3*k);
                bphi = bphi + (dR_dphi * Cs((3*ind_cont(j)-2):3*ind_cont(j),(3*k-2):3*k) * R' + R * Cs((3*ind_cont(j)-2):3*ind_cont(j),(3*k-2):3*k) * dR_dphi') * FSurf3((3*k-2):3*k);
                bpsi = bpsi + (dR_dpsi * Cs((3*ind_cont(j)-2):3*ind_cont(j),(3*k-2):3*k) * R' + R * Cs((3*ind_cont(j)-2):3*ind_cont(j),(3*k-2):3*k) * dR_dpsi') * FSurf3((3*k-2):3*k);
                W_Fs = W_Fs + Ws((3*ind_cont(j)-2):3*ind_cont(j),(3*k-2):3*k) * dF_dPfree_surf((3*k-2):3*k,(3*g-2):3*g); 
                RdCRt1 = RdCRt1 + R * der_C_s{(3*g)-2}((3*ind_cont(j)-2):3*ind_cont(j),(3*k-2):3*k) * R' * FSurf3((3*k-2):3*k);
                RdCRt2 = RdCRt2 + R * der_C_s{(3*g)-1}((3*ind_cont(j)-2):3*ind_cont(j),(3*k-2):3*k) * R' * FSurf3((3*k-2):3*k);
                RdCRn = RdCRn + R * der_C_s{3*g}((3*ind_cont(j)-2):3*ind_cont(j),(3*k-2):3*k) * R' * FSurf3((3*k-2):3*k);
            end
            A = zeros(3,3);
            if ~isempty(find(g==ind_cont,1))
                if (j==contj && flag==0)
                   A = R;
                   contj = contj + 1;
                   flag = 1;
                end
            end
            dPc_dPfrees((3*j-2):3*j,(3*g-2):3*g) = dOl_dPfree_surf(:,(3*g-2):3*g) + [btheta bphi bpsi] * dY_dPfree_surf(:,(3*g-2):3*g) + A + W_Fs + [RdCRt1 RdCRt2 RdCRn];
        end
        flag=0;
        
    end
    der = dPc_dPfrees * dPFreeSurf3_dp;
%     [der(:,9) ((a.P_mc3-P_mc3)./1e-09)]
    %dP_dPfree = H_tilde * Ccc * RBlock_c' + RBlock_c * der_C_cc * RBlock_c' * Fc_mc3 + Wc * dF_dP + dOl_dP + H_tilde * P_mc3 + RBlock_c;
    
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %    dPInt_dP_i              %
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     dPc_dPInt = zeros(3*length(ind_cont),3*sole.nInt);
%     [dR_dtheta,dR_dphi,dR_dpsi] = dRot(angleact);
%     contj = 1;
%     flag=0;
%     for g=1:sole.nInt
%         for j=1:length(ind_cont)
% %             btheta = dR_dtheta*Pfree_c3((3*j-2):3*j);
% %             bphi = dR_dphi*Pfree_c3((3*j-2):3*j);
% %             bpsi = dR_dpsi*Pfree_c3((3*j-2):3*j);
%             W_Fs = zeros(3,3);
%             RdCRn = zeros(3,1);
%             RdCRt1 = zeros(3,1);
%             RdCRt2 = zeros(3,1);
%             for k=1:sole.nFreeSurf
% %                 btheta = btheta + (dR_dtheta * sole.WC_ii_is((3*ind_cont(j)-2):3*ind_cont(j),(3*k-2):3*k) * R' + R * sole.WC_ii_is((3*ind_cont(j)-2):3*ind_cont(j),(3*k-2):3*k) * dR_dtheta') * FSurf3((3*k-2):3*k);
% %                 bphi = bphi + (dR_dphi * Cs((3*ind_cont(j)-2):3*ind_cont(j),(3*k-2):3*k) * R' + R * Cs((3*ind_cont(j)-2):3*ind_cont(j),(3*k-2):3*k) * dR_dphi') * FSurf3((3*k-2):3*k);
% %                 bpsi = bpsi + (dR_dpsi * Cs((3*ind_cont(j)-2):3*ind_cont(j),(3*k-2):3*k) * R' + R * Cs((3*ind_cont(j)-2):3*ind_cont(j),(3*k-2):3*k) * dR_dpsi') * FSurf3((3*k-2):3*k);
% %                 W_Fs = W_Fs + Ws((3*ind_cont(j)-2):3*ind_cont(j),(3*k-2):3*k) * dF_dPfree_surf((3*k-2):3*k,(3*g-2):3*g); 
%                 RdCRt1 = RdCRt1 + R * sole.der_m_invKii_Kis{(3*g)-2}((3*ind_cont(j)-2):3*ind_cont(j),(3*k-2):3*k) * R' * FSurf3((3*k-2):3*k);
%                 RdCRt2 = RdCRt2 + R * sole.der_m_invKii_Kis{(3*g)-1}((3*ind_cont(j)-2):3*ind_cont(j),(3*k-2):3*k) * R' * FSurf3((3*k-2):3*k);
%                 RdCRn = RdCRn + R * sole.der_m_invKii_Kis{3*g}((3*ind_cont(j)-2):3*ind_cont(j),(3*k-2):3*k) * R' * FSurf3((3*k-2):3*k);
%             end
%             A = zeros(3,3);
%             if ~isempty(find(g==ind_cont,1))
%                 if (j==contj && flag==0)
%                    A = R;
%                    contj = contj + 1;
%                    flag = 1;
%                 end
%             end
%             dPc_dPfrees((3*j-2):3*j,(3*g-2):3*g) = A + [RdCRt1 RdCRt2 RdCRn];
%         end
%         flag=0;
%         
%     end
%     der = dPc_dPfrees * dPFreeSurf3_dp;

%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%     %    dPsurf_dPfrees          %
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     dPsurf_dPfrees = zeros(3*sole.nFreeSurf,3*sole.nFreeSurf);
%     [dR_dtheta,dR_dphi,dR_dpsi] = dRot(angleact);
%     for g=1:sole.nFreeSurf
%         for j=1:sole.nFreeSurf
%             btheta = dR_dtheta*Pfree_s3((3*j-2):3*j);
%             bphi = dR_dphi*Pfree_s3((3*j-2):3*j);
%             bpsi = dR_dpsi*Pfree_s3((3*j-2):3*j);
%             W_Fs = zeros(3,3);
%             RdCRn = zeros(3,1);
%             RdCRt1 = zeros(3,1);
%             RdCRt2 = zeros(3,1);
%             for k=1:sole.nFreeSurf
%                 btheta = btheta + (dR_dtheta * Cs((3*j-2):3*j,(3*k-2):3*k) * R' + R * Cs((3*j-2):3*j,(3*k-2):3*k) * dR_dtheta') * FSurf3((3*k-2):3*k);
%                 bphi = bphi + (dR_dphi * Cs((3*j-2):3*j,(3*k-2):3*k) * R' + R * Cs((3*j-2):3*j,(3*k-2):3*k) * dR_dphi') * FSurf3((3*k-2):3*k);
%                 bpsi = bpsi + (dR_dpsi * Cs((3*j-2):3*j,(3*k-2):3*k) * R' + R * Cs((3*j-2):3*j,(3*k-2):3*k) * dR_dpsi') * FSurf3((3*k-2):3*k);
%                 W_Fs = W_Fs + Ws((3*j-2):3*j,(3*k-2):3*k) * dF_dPfree_surf((3*k-2):3*k,(3*g-2):3*g); 
%                 RdCRt1 = RdCRt1 + R * der_C_s{(3*g)-2}((3*j-2):3*j,(3*k-2):3*k) * R' * FSurf3((3*k-2):3*k);
%                 RdCRt2 = RdCRt2 + R * der_C_s{(3*g)-1}((3*j-2):3*j,(3*k-2):3*k) * R' * FSurf3((3*k-2):3*k);
%                 RdCRn = RdCRn + R * der_C_s{3*g}((3*j-2):3*j,(3*k-2):3*k) * R' * FSurf3((3*k-2):3*k);
%             end
%             A=zeros(3,3);
%             if j==g
%                A=R;
%             end
%             dPsurf_dPfrees((3*j-2):3*j,(3*g-2):3*g) = dOl_dPfree_surf(:,(3*g-2):3*g) + [btheta bphi bpsi] * dY_dPfree_surf(:,(3*g-2):3*g) + A + W_Fs + [RdCRt1 RdCRt2 RdCRn];
%         end
%     end
    %dP_dPfree = H_tilde * Ccc * RBlock_c' + RBlock_c * der_C_cc * RBlock_c' * Fc_mc3 + Wc * dF_dP + dOl_dP + H_tilde * P_mc3 + RBlock_c;
    
    %%% Find the displacement of internal nodes and position of all nodes
    FcSurf3 = zeros(D*sole.nFreeSurf,1);
    FcSurf3(ind_cont3) = Fc_mc3;
    FcFreeSurf = reshape(FcSurf3,D,m)';
    Fcloc = R'*FcFreeSurf';
    Fcloc3 = reshape(Fcloc,D*m,1);
    dPlocSurf3 = sole.Cs * Fcloc3;
    % Displacement of all nodes
    dPInt3 = -(sole.m_invKii_Kis) * dPlocSurf3;
    dPloc3 = zeros(D*sole.nTot,1);
    dPloc3(sole.nodesFreeSurf3) = dPlocSurf3;
    dPloc3(sole.nodesInt3) = dPInt3;
    dPloc = reshape(dPloc3,D,sole.nTot)';
    dPabstot = R*dPloc';
    displacement = repmat(displ', sole.nTot, 1);
    % PosRot = R*P0 + center;
    PosRot = R*sole.coor';
    Pg = PosRot' + dPabstot' + displacement;
    %Pg(sole.nodesFreeSurf,:)'=PSurf;
    Fc = FcFreeSurf(ind_cont,:);
    Pc = Pg(sole.nodesFreeSurf(ind_cont),:);
    
    stressVM0 = zeros(sole.nTot,1); % VonMises' Stress    
    plotsole(2,sole.elements_surf,Pg,stressVM0,Pc,Fc,Z,Ftot,-127.5,30);

end
P_mc3_surf0 = [0.000000000000000;
   0.064969879067856;
  -0.000000000000000;
   0.000000000000000;
   0.012993975813440;
   0.000000000000000;
   0.000000000000000;
   0.038981927440646;
  -0.000000000000000;
   0.028741597655377;
   0.064991035607552;
  -0.000000000000000;
   0.057483195310722;
   0.065012192147248;
  -0.000000000000000;
   0.086542434327247;
   0.064678457534719;
  -0.000000000000000;
   0.020196927969467;
   0.044984370312278;
  -0.000000000000000;
   0.015184423223862;
   0.025160715975528;
   0.000000000000000;
   0.041622393494572;
   0.049145268495376;
  -0.000000000000000];
P_mc3_surf1 = [-0.000000000000000;
   0.064969879067856;
  -0.000000000000000;
  -0.000000000000000;
   0.012993975813440;
  -0.000000000000000;
  -0.000000000000000;
   0.038981927440646;
   0.000000000000000;
   0.028741597655377;
   0.064991035607552;
   0.000000000000000;
   0.057483195310722;
   0.065012192147248;
   0.000000000000000;
   0.086542434327247;
   0.064678457534719;
   0.000000000000000;
   0.020196927969467;
   0.044984370312278;
  -0.000000000000000;
   0.015184423223861;
   0.025160715975529;
  -0.000000000000000;
   0.041622393494572;
   0.049145268495376;
   0.000000000000000];

P_mc3_1 = [0.000000000999708;
   0.064969879067757;
  -0.000000000000000;
   0.000000000000000;
   0.012993975813273;
   0.000000000000000;
   0.000000000000000;
   0.038981927440145;
  -0.000000000000000;
   0.028741597655180;
   0.064991035607259;
   0.000000000000000;
   0.057483195310327;
   0.065012192147497;
   0.000000000000000;
   0.086542434333162;
   0.064678457521922;
  -0.000000000000000;
   0.020196927969328;
   0.044984370312080;
  -0.000000000000000;
   0.015184423224795;
   0.025160715973504;
   0.000000000000000;
   0.041622393497287;
   0.049145268492831;
  -0.000000000000000];
P_mc3_0 = [
   0.000000000000000;
   0.064969879067856;
  -0.000000000000000;
   0.000000000000000;
   0.012993975813440;
   0.000000000000000;
   0.000000000000000;
   0.038981927440646;
  -0.000000000000000;
   0.028741597655377;
   0.064991035607552;
  -0.000000000000000;
   0.057483195310722;
   0.065012192147248;
  -0.000000000000000;
   0.086542434327247;
   0.064678457534719;
  -0.000000000000000;
   0.020196927969467;
   0.044984370312278;
  -0.000000000000000;
   0.015184423223862;
   0.025160715975528;
   0.000000000000000;
   0.041622393494572;
   0.049145268495376;
  -0.000000000000000];

Fc_mc3_1 = [5.277001438200586;
  -5.770497604403332;
  16.330217824383812;
   0.934631089615455;
   1.763854775915847;
   5.256199415130856;
   5.848392758215681;
   0.683924628450662;
  13.881977712784307;
  -1.121459834020889;
  -7.433786599032373;
  17.782075330536475;
  -3.475016546138555;
  -1.090234671500249;
  10.709890159139519;
  -0.008187474545790;
   0.009147619555551;
   0.015345687656220;
  -0.923576797396979;
   1.745719685824719;
  17.761583849541843;
  -4.602214497058615;
   7.723057614891084;
  11.237909975826767;
  -6.528386136870724;
   5.052287550298045;
  10.318781045000428];

Fc_mc3 = [5.277001414567717;
  -5.770497472165673;
  16.330217598756825;
   0.934631080797675;
   1.763854772341003;
   5.256199425673882;
   5.848392684666790;
   0.683924693684816;
  13.881977452984538;
  -1.121459803377909;
  -7.433786677097638;
  17.782075491161372;
  -3.475016479098905;
  -1.090234739129581;
  10.709890125211242;
  -0.008187531079164;
   0.009147682555784;
   0.015345793464419;
  -0.923576660231524;
   1.745719572581303;
  17.761584022229329;
  -4.602214541296511;
   7.723057653971813;
  11.237910046098815;
  -6.528386164947742;
   5.052287513258059;
  10.318781044419106];
P_mc3_p1=[  -0.000000000000179;
   0.064969879067963;
   0.000000000000000;
   0.000000000000000;
   0.012993975813441;
   0.000000000000000;
  -0.000000000001313;
   0.038981927441087;
   0.000000000000000;
   0.028741597655382;
   0.064991035607551;
   0.000000000000000;
   0.057483195310731;
   0.065012192147239;
   0.000000000000000;
   0.086542434327065;
   0.064678457535084;
  -0.000000000000000;
   0.020196927969231;
   0.044984370312389;
   0.000000000000000;
   0.015184423223913;
   0.025160715975423;
   0.000000000000000;
   0.041622393494545;
   0.049145268495396;
   0.000000000000000];


P_mc3_p0 = [-0.000000000000000;
   0.064969879067856;
  -0.000000000000000;
   0.000000000000000;
   0.012993975813440;
   0.000000000000000;
   0.000000000000000;
   0.038981927440646;
  -0.000000000000000;
   0.028741597655377;
   0.064991035607552;
   0.000000000000000;
   0.057483195310722;
   0.065012192147248;
   0.000000000000000;
   0.086542434327247;
   0.064678457534719;
  -0.000000000000000;
   0.020196927969467;
   0.044984370312278;
                   0;
   0.015184423223862;
   0.025160715975528;
   0.000000000000000;
   0.041622393494572;
   0.049145268495376;
  -0.000000000000000];

G_Fc = (Fc_mc3_1-Fc_mc3)/(1e-09);
E_rel=(dF_ddelta(1:3*length(ind_cont),1) - G_Fc)./G_Fc;
[dF_ddelta(1:3*length(ind_cont),1) G_Fc]
displ1 = 1.0e-03 * [0.249008261247135;-0.275332293301581;-0.141828523015672];
displ0 = 1.0e-03 * [0.249008264182643;-0.275332292903205;-0.141828530268763];
G1 = (displ1-displ0)./1e-09;
angleact1 = [-0.030444531958799;-0.024177249920720;-0.001910632277161];
angleact0 = [-0.030444531536849;-0.024177249636684;-0.001910632163002];
G1 = (angleact1-angleact0)./1e-09;
E_rel = (G - G1)./abs(G1);

%%%%% 1 e-06
displ1 = 1.0e-03 * [0.249005328570552;-0.275332691251498;-0.141821277063753];
displ0 = 1.0e-03 * [0.249008264182643;-0.275332292903205;-0.141828530268763]; 