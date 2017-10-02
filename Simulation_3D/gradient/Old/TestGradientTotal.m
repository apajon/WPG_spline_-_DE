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
pname = 'input/semelle1 L=0.23, l=0.13, e=0.03 m new centre/';
%%% foot size %%%
l=0.13;
L=0.23;
e=0.03;
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
[sole.coor,dPtot3_dp] = deformation_moveDiri(sole,p_ini_v,spl,move_dirichlet);
dPFreeSurf3_dp = dPtot3_dp(sole.nodesFreeSurf3,:);

st1 = 1e-09;
% sole.coor(2,1) = sole.coor(2,1) + st1;
% sole.coor(131,1) = sole.coor(131,1) + st1;


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
%     Asurf = gradient_simu_A_surf(sole,P_mc,Fc_mc,ind_cont,Ws,ind_slip,PabsOld_mc,fric,contAngle,displ_first,psi_first,angleact);
% 
    Ccc = Cs(ind_cont3,ind_cont3);
    B = gradient_simu_B_contact(angleact,Pfree_c3,Fc_mc3,ind_cont,Ccc,ind_slip,psi_first,contAngle);
    Pfree_s = Pfree;
    Pfree_s3 = reshape(Pfree_s,3*sole.nFreeSurf,1);
    FSurf3 = zeros(D*sole.nFreeSurf,1);
    FSurf3(ind_cont3) = Fc_mc3;
%     FSurf = reshape(FcSurf3,D,m)';
%     Bsurf = gradient_simu_B_surf(sole,angleact,Pfree_s3,FSurf3,Cs,ind_slip,psi_first,contAngle);
    A_alg = A_out(1:(D*length(ind_cont)+D*length(ind_slip))*(D*length(ind_cont)+D*length(ind_slip)));
    A_alg = reshape(A_alg,D*length(ind_cont)+D*length(ind_slip),D*length(ind_cont)+D*length(ind_slip));
    B_alg = B_out(1:6*(D*length(ind_cont)+D*length(ind_slip)));
    B_alg = reshape(B_alg,D*length(ind_cont)+D*length(ind_slip),6);
%     G_alg = G_out(1:6*(D*length(ind_cont)+D*length(ind_slip)));
%     G_alg = reshape(G_alg,6,D*length(ind_cont)+D*length(ind_slip));
    J2_alg = J2(1:(6*6));
    J2_alg = reshape(J2_alg,6,6);
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
    Theta_mat_surf = gradient_simu_Theta_surf(sole,der_C_s,ind_slip,FSurf3,RBlock_s,Rini_block_s);
    E = gradient_simu_E_contact(ind_cont,ind_slip,Fc_mc,P_mc,Zdes);
%     Esurf = gradient_simu_E_surf(sole,ind_cont,ind_slip,Fc_mc,P_mc,Zdes);
    N_ini = gradient_simu_N_ini_contact(ind_cont,Rini,Pfree_c3,Fc_mc3,angleact,psi_first,contAngle);
%     N_ini_surf = gradient_simu_N_ini_surf(ind_cont,Rini,Pfree_c3,Fc_mc3,angleact,psi_first,contAngle);
    Y = gradient_simu_Y_contact(ind_cont,Rini,Fc_mc3);
%     Nini = gradient_simu_Nini_contact(ind_cont,Rini,P_mc3,Fc_mc3,angleact,psi_first,contAngle);
    %G1 = inv(A) * Theta_mat;
%     Grad1 = inv((E / A * B) + N_ini)*((-E / A *Theta_mat) + Y);
%     %Grad1 = inv((E / A *B) + Nini)*((-E / A *Theta_mat));
%     Gc = Grad1(1:3,1);
    inv_A = inv(A);
    A_1B = inv_A * B;
    inv_A_theta = inv_A * Theta_mat;
    dOl_dY = inv((E * A_1B) + N_ini)*((-E * inv_A_theta) + Y);
    dOl_dPfree_c = dOl_dY(1:3,:);
    dOl_dPfree_surf = zeros(3,3*sole.nFreeSurf);
    dOl_dPfree_surf(:,ind_cont3) = dOl_dPfree_c;
    dOl_dPfree = zeros(3,3*sole.nFree);
    dOl_dPfree(:,sole.nodesFreeSurf3) = dOl_dPfree_surf;
    dOlsurf3_dp = dOl_dPfree_surf * dPFreeSurf3_dp;
    dY_dPfree_c = dOl_dY(4:6,:);
    dY_dPfree_surf = zeros(3,3*sole.nFreeSurf);
    dY_dPfree_surf(:,ind_cont3) = dY_dPfree_c;
    dY_dPfree = zeros(3,3*sole.nFree);
    dY_dPfree(:,sole.nodesFreeSurf3) = dY_dPfree_surf;
    dYsurf3_dp = dY_dPfree_surf * dPFreeSurf3_dp;
    %Gc = dOl_dY(1:3,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %    dF_dPfree and ddelta_dPfree     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dF_ddeltat = A_1B * dOl_dY + inv_A_theta;
    dF_dPfree_c = dF_ddeltat(1:3*length(ind_cont),:);
    dF_dPfree_surf = zeros(3*sole.nFreeSurf,3*sole.nFreeSurf);
    dFsurf3_dp = dF_dPfree_surf * dPFreeSurf3_dp;
    for j=1:3*length(ind_cont)
        dF_dPfree_surf(ind_cont3,ind_cont3(j)) = dF_dPfree_c(:,j);
    end
    dF_dPfree = zeros(3*sole.nFree,3*sole.nFree);
    for j=1:3*sole.nFreeSurf
        dF_dPfree(sole.nodesFreeSurf3,sole.nodesFreeSurf3(j)) = dF_dPfree_surf(:,j);
    end
    ddeltat_dPfree = dF_ddeltat(3*length(ind_cont)+1:end,:);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %    dPInt_dP_i              %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     dPInt_dPfree_is = zeros(3*sole.nFree,3*sole.nInt);
%     [dR_dtheta,dR_dphi,dR_dpsi] = dRot(angleact);
%     for g=1:sole.nFree
%         for j=1:sole.nInt
%             btheta = dR_dtheta*Pfree_s3((3*j-2):3*j);
%             bphi = dR_dphi*Pfree_s3((3*j-2):3*j);
%             bpsi = dR_dpsi*Pfree_s3((3*j-2):3*j);
%             RdCRn = zeros(3,1);
%             RdCRt1 = zeros(3,1);
%             RdCRt2 = zeros(3,1);
%             W_Fs = zeros(3,3);
%             otheta = zeros(3,1);
%             ophi = zeros(3,1);
%             opsi = zeros(3,1);
%             for k=1:sole.nFreeSurf
%                 otheta = otheta + (dR_dtheta * sole.WC_ii_is((3*j-2):3*j,(3*k-2):3*k) * R' + R * sole.WC_ii_is((3*j-2):3*j,(3*k-2):3*k) * dR_dtheta') * FSurf3((3*k-2):3*k);
%                 ophi = ophi + (dR_dphi * sole.WC_ii_is((3*j-2):3*j,(3*k-2):3*k) * R' + R * sole.WC_ii_is((3*j-2):3*j,(3*k-2):3*k) * dR_dphi') * FSurf3((3*k-2):3*k);
%                 opsi = opsi + (dR_dpsi * sole.WC_ii_is((3*j-2):3*j,(3*k-2):3*k) * R' + R * sole.WC_ii_is((3*j-2):3*j,(3*k-2):3*k) * dR_dpsi') * FSurf3((3*k-2):3*k);
%                 W_Fs = W_Fs + R * sole.WC_ii_is((3*j-2):3*j,(3*k-2):3*k) * R' * dF_dPfree((3*k-2):3*k,(3*g-2):3*g);
%                 RdCRt1 = RdCRt1 + R * sole.der_m_invKii_Kis{(3*g)-2}((3*j-2):3*j,(3*k-2):3*k) * R' * FSurf3((3*k-2):3*k);
%                 RdCRt2 = RdCRt2 + R * sole.der_m_invKii_Kis{(3*g)-1}((3*j-2):3*j,(3*k-2):3*k) * R' * FSurf3((3*k-2):3*k);
%                 RdCRn = RdCRn + R * sole.der_m_invKii_Kis{3*g}((3*j-2):3*j,(3*k-2):3*k) * R' * FSurf3((3*k-2):3*k);
%             end
%             A=zeros(3,3);
% %             if j==g
% %                A=R;
% %             end
%             dPInt_dPfree_is((3*j-2):3*j,(3*k-2):3*k) = [otheta ophi opsi] * dY_dPfree_surf(:,(3*j-2):3*j) + [RdCRt1 RdCRt2 RdCRn] + W_Fs;
%             %dPInt_dPfree_is((3*g-2):3*g,(3*j-2):3*j) = [RdCRt1 RdCRt2 RdCRn];
%         end
%     end
    
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
%     %dP_dPfree = H_tilde * Ccc * RBlock_c' + RBlock_c * der_C_cc * RBlock_c' * Fc_mc3 + Wc * dF_dP + dOl_dP + H_tilde * P_mc3 + RBlock_c;


%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%     %    dPc_dPfrees             %
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     dPsurf_dPfrees = zeros(3*sole.nFreeSurf,3*length(ind_cont));
%     [dR_dtheta,dR_dphi,dR_dpsi] = dRot(angleact);
%     for g=1:length(ind_cont)
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
%                 W_Fs = W_Fs + Ws((3*j-2):3*j,(3*k-2):3*k) * dF_dPfree_surf((3*k-2):3*k,(3*sole.nodesFreeSurf3(g)-2):3*sole.nodesFreeSurf3(g)); 
%                 RdCRt1 = RdCRt1 + R * der_C_s{(3*sole.nodesFreeSurf3(g))-2}((3*j-2):3*j,(3*k-2):3*k) * R' * FSurf3((3*k-2):3*k);
%                 RdCRt2 = RdCRt2 + R * der_C_s{(3*sole.nodesFreeSurf3(g))-1}((3*j-2):3*j,(3*k-2):3*k) * R' * FSurf3((3*k-2):3*k);
%                 RdCRn = RdCRn + R * der_C_s{3*sole.nodesFreeSurf3(g)}((3*j-2):3*j,(3*k-2):3*k) * R' * FSurf3((3*k-2):3*k);
%             end
%             A=zeros(3,3);
%             if j==sole.nodesFreeSurf3(g)
%                A=R;
%             end
%             dPsurf_dPfrees((3*j-2):3*j,(3*g-2):3*g) = dOl_dPfree_surf(:,(3*g-2):3*g) + [btheta bphi bpsi] * dY_dPfree_surf(:,(3*g-2):3*g) + A + W_Fs + [RdCRt1 RdCRt2 RdCRn];
%         end
%     end
    %dP_dPfree = H_tilde * Ccc * RBlock_c' + RBlock_c * der_C_cc * RBlock_c' * Fc_mc3 + Wc * dF_dP + dOl_dP + H_tilde * P_mc3 + RBlock_c;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %    dPc_dPfrees             %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     dPc_dPfree = zeros(3*sole.nTot,3*length(ind_cont));
%     cont = 1;    
%     for j=1:3*length(ind_cont)    
%             contint = 1;
%             for k=1:(3*sole.nTot)
%                 if contint <= 3*length(ind_cont)
%                     if k==sole.nodesFreeSurf3(ind_cont3(contint))
%                         dPc_dPfree(k,ind_cont3(j)) = dPc_dPfreec(contint,j);
%                         contint = contint + 1;
%                     end
%                 end
%             end
%     end    
    
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
%     Psint = [Pg(sole.nodesFreeSurf,:); Pg(sole.nodesInt,:)];
%     Psint3 = reshape(Psint',3*(sole.nFreeSurf+sole.nInt),1);
%     save 'Pint3.mat' Psint3
    
    stressVM0 = zeros(sole.nTot,1); % VonMises' Stress    
    plotsole(2,sole.elements_surf,Pg,stressVM0,Pc,Fc,Z,Ftot,-127.5,30);

end

Pint_0 = [-0.004772835698615;
  -0.029611515116526;
   0.036496548052675;
  -0.003719222530341;
   0.011983393983610;
   0.033091653028863;
  -0.004110702277668;
   0.040923043353869;
   0.030169079768867;
  -0.009248338495276;
  -0.007766371989636;
   0.044672121742295;
   0.001409086783416;
  -0.008938568080645;
   0.024154997386469;
   0.009931475761847;
   0.030162704047845;
   0.041714181726227];

Pint_1 = [   0.029639115747695;
  -0.058441693034852;
   0.062922219157346;
   0.018150959870111;
  -0.063002462130303;
   0.004159017174040;
   0.032304359172698;
   0.071377371279591;
   0.055002484769576;
   0.026195173449623;
   0.064640590962397;
  -0.000000000000000;
   0.025814529910867;
  -0.059972928548122;
   0.043335863857261;
   0.021971317310846;
  -0.061496450390414;
   0.023747325067460;
  -0.001499546414464;
  -0.062210069343564;
   0.007964500151788;
  -0.021236777985021;
  -0.061549661740905;
   0.011872595566548;
  -0.009431129345423;
  -0.056505343419225;
   0.070501968055748;
   0.010046997132188;
  -0.057550615062542;
   0.066687819263374;
   0.029578261067549;
   0.069872222559673;
   0.035579811630067;
   0.027086823500616;
   0.068351166256511;
   0.017045831122287;
   0.004858094545349;
   0.066594612570891;
  -0.000000000000000;
  -0.016896684542811;
   0.067908445612192;
   0.000343699659645;
  -0.006869900223079;
   0.073169270181481;
   0.059253451974144;
   0.012543348448685;
   0.072174498848519;
   0.057322519802567;
   0.030218846008515;
  -0.036861850175480;
   0.061200314056290;
   0.030674711325045;
  -0.015358037030185;
   0.059666597671457;
   0.031023269886070;
   0.006153276279958;
   0.058453286224461;
   0.031279859080315;
   0.027764572585311;
   0.057604329403077;
   0.031751469699092;
   0.049451988146921;
   0.056753731752594;
   0.018797782174023;
  -0.041369844219499;
   0.002369591630797;
   0.019667569804126;
  -0.019572917967455;
   0.000686126469472;
   0.021062487126294;
   0.002991315872569;
  -0.000000000000000;
   0.022360137415689;
   0.025313372772274;
  -0.000000000000000;
   0.023734497044114;
   0.047145402630235;
  -0.000000000000000;
  -0.005549109084846;
  -0.059478737276494;
   0.039182603518299;
  -0.015288131984925;
  -0.057801146811598;
   0.053820826266380;
   0.004316599657464;
  -0.060958559321664;
   0.024592893616803;
  -0.020304779574183;
  -0.060056574551090;
   0.029427146224796;
   0.009123222438417;
  -0.059047093964531;
   0.049068189222525;
  -0.024783032563262;
  -0.058601935886979;
   0.043060132323825;
   0.013732325916823;
  -0.060300211549309;
   0.035497964899662;
  -0.009364999099133;
  -0.061048220407396;
   0.020014573879633;
  -0.001709059166300;
  -0.057876676282008;
   0.058437662429679;
  -0.002190393061626;
   0.070186674878667;
   0.028757579977418;
  -0.012136706427365;
   0.071752244228515;
   0.042088045182824;
   0.008539220601828;
   0.068780388619129;
   0.016099761570453;
  -0.016977276489472;
   0.069389103211129;
   0.017718380798339;
   0.012343547156033;
   0.070629171380201;
   0.039817728970095;
  -0.021626065232758;
   0.070809539908738;
   0.030668314574843;
   0.017727968922388;
   0.069550444662622;
   0.027248433808987;
  -0.005246883774423;
   0.068468621721865;
   0.009912954562905;
   0.001240098691442;
   0.071794954110701;
   0.048027888455704;
  -0.010493833359928;
  -0.032655202380543;
   0.007285414625713;
  -0.007321834107095;
   0.037108898400423;
   0.001762724308969;
  -0.009296758395919;
  -0.002152984801362;
   0.004722840269521;
  -0.016667173024894;
   0.018929772047433;
   0.004291004283634;
   0.003478609400678;
   0.015994521827000;
   0.001583103957823;
  -0.021014117832044;
  -0.020207206896986;
   0.008193968093194;
   0.000849505546850;
  -0.021773304579215;
   0.004293224452085;
  -0.023698672053186;
  -0.044189792158313;
   0.010877250927528;
   0.001851508086094;
  -0.045031459673890;
   0.005889839686772;
  -0.021929018696680;
   0.051958127133852;
   0.002145517490246;
   0.006487677987759;
   0.050553974601217;
   0.000216243696166;
  -0.024754705833161;
   0.035974102841289;
   0.003807610939489;
   0.009552848896519;
   0.034808232705876;
   0.000298205068145;
  -0.023598937305823;
   0.004159159415447;
   0.006597490711679;
   0.006532444151035;
  -0.003608805754326;
   0.002175212724310;
  -0.011042208305769;
  -0.049122518505593;
   0.008765427141785;
  -0.007056994779317;
   0.057085791266032;
   0.000591861894889;
   0.000991273200936;
  -0.027894596768783;
   0.066009110762114;
   0.001680878394409;
   0.041770022897406;
   0.060810138683715;
   0.001454993188321;
   0.002282457221226;
   0.063595375141337;
   0.009648369428046;
   0.022709692831023;
   0.060964346073523;
  -0.009871594543772;
   0.021257207925825;
   0.063955108737117;
   0.011927315934944;
  -0.016381394239970;
   0.063023546827194;
  -0.009560880488942;
  -0.016153683827235;
   0.066907544246263;
   0.013370539883721;
  -0.040408534824699;
   0.064603602742399;
  -0.011829694190511;
  -0.039058934406071;
   0.069467637284650;
   0.016407705204736;
   0.055764828499522;
   0.058129543657534;
  -0.010636043112984;
   0.056314561082096;
   0.060943050620398;
   0.017946025653011;
   0.039381886581298;
   0.058881236483991;
  -0.014392908957438;
   0.040466408498453;
   0.062775763042096;
   0.015655301975684;
   0.007498804250885;
   0.060894813716282;
  -0.013918819268154;
   0.002142385141340;
   0.066087475973316;
   0.000644879640083;
  -0.044293062864905;
   0.067435931981053;
   0.002585433745794;
   0.062007051498824;
   0.059140360152134;
   0.024690275246728;
  -0.031538451127683;
   0.031199439806220;
   0.027160142714470;
   0.038101435176533;
   0.028373166539517;
   0.025663417258199;
  -0.001226064788523;
   0.029254055485442;
   0.024971706156258;
   0.019246348202583;
   0.020878184966116;
   0.028213133760231;
   0.017983305631538;
   0.040288882783179;
   0.023003630190050;
  -0.020240694364052;
   0.019466763971123;
   0.027065108556877;
  -0.019570477238187;
   0.041243090161745;
   0.021879281842014;
  -0.044432722857767;
   0.019433334442891;
   0.026842099357622;
  -0.042456909034152;
   0.044778607154206;
   0.025527883559181;
   0.051707689610276;
   0.014247534874231;
   0.029508055289426;
   0.053059345348693;
   0.040318280945003;
   0.024551170679496;
   0.035857626053833;
   0.012313217614022;
   0.029594495646684;
   0.037309185435499;
   0.044606934933662;
   0.023421828535691;
   0.003888572728606;
   0.014923856235124;
   0.028414240566910;
  -0.001129800988800;
   0.044910936461255;
   0.024240479658145;
  -0.047981440303141;
   0.032484148275157;
   0.027816700693277;
   0.058473333016062;
   0.027096719039507];

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
P_mc3_0 = [-0.000000000000000;
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




P_mc3_int1 =[-0.000000000000000;
   0.064969409453513;
   0.000000000000000;
   0.000036034759253;
   0.007140849975361;
  -0.000000000000000;
  -0.000000000000000;
   0.021656469817727;
   0.000000000000000;
  -0.000000000000000;
   0.036094116362986;
   0.000000000000000;
  -0.000000000000000;
   0.050531762908245;
   0.000000000000000;
   0.014370376136650;
   0.064980593200307;
   0.000000000000000;
   0.028740752273307;
   0.064991776947101;
   0.000000000000000;
   0.043111128409945;
   0.065002960693896;
  -0.000000000000000;
   0.057481504546582;
   0.065014144440690;
   0.000000000000000;
   0.071897752387106;
   0.065014619321365;
   0.000000000000000;
   0.033772337109929;
   0.038651452391521;
  -0.000000000000000;
   0.053248316770119;
   0.044858697717188;
   0.000000000000000;
   0.019399004183838;
   0.024702741416337;
  -0.000000000000000;
   0.020357476117264;
   0.044533647104384;
   0.000000000000000;
   0.037295848404654;
   0.053960791499794;
   0.000000000000000;
   0.011528675868127;
   0.036858357673683;
  -0.000000000000000;
   0.050108153172051;
   0.055887249683516;
   0.000000000000000;
   0.011955032841029;
   0.052770105246252;
   0.000000000000000;
   0.023405226508995;
   0.034497094149880;
  -0.000000000000000;
   0.025498175428667;
   0.054722227583377;
   0.000000000000000;
   0.043562717828692;
   0.044551995349102;
   0.000000000000000];


P_mc3_int0 =[-0.000000000000000;
   0.064969409453513;
   0.000000000000000;
   0.000036034759253;
   0.007140849975361;
  -0.000000000000000;
  -0.000000000000000;
   0.021656469817727;
   0.000000000000000;
  -0.000000000000000;
   0.036094116362986;
   0.000000000000000;
  -0.000000000000000;
   0.050531762908245;
   0.000000000000000;
   0.014370376136650;
   0.064980593200307;
   0.000000000000000;
   0.028740752273307;
   0.064991776947101;
   0.000000000000000;
   0.043111128409945;
   0.065002960693896;
  -0.000000000000000;
   0.057481504546582;
   0.065014144440690;
   0.000000000000000;
   0.071897752387106;
   0.065014619321365;
   0.000000000000000;
   0.033772337109929;
   0.038651452391521;
  -0.000000000000000;
   0.053248316770119;
   0.044858697717188;
   0.000000000000000;
   0.019399004183838;
   0.024702741416337;
  -0.000000000000000;
   0.020357476117264;
   0.044533647104384;
   0.000000000000000;
   0.037295848404654;
   0.053960791499794;
   0.000000000000000;
   0.011528675868127;
   0.036858357673683;
  -0.000000000000000;
   0.050108153172051;
   0.055887249683516;
   0.000000000000000;
   0.011955032841029;
   0.052770105246252;
   0.000000000000000;
   0.023405226508995;
   0.034497094149880;
  -0.000000000000000;
   0.025498175428667;
   0.054722227583377;
   0.000000000000000;
   0.043562717828692;
   0.044551995349102;
   0.000000000000000];


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