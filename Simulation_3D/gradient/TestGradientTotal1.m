format long
clc
% clear all
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
sole = soleFEM_newStiff(pname,fname,l,L,e);
soleini = soleFEM_newStiff(pname,fname,l,L,e);
coorini = sole.coor;
Young = 1000000;
Poisson = 0.3;
sole.setMaterial(Young,Poisson);
Fdes = 1.0e+02 * [-0.052709310000000;0.030650450000000;1.279851840000000];
Zdes = [0.020625645823457;0.049298000000611];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           B-spline                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
spline_res = -pi/2:pi/10:pi/2;
l_spli = length(spline_res);
spl = SplineClass(sole,spline_res);
%splini = SplineClass(soleini,spline_res);

move_dirichlet = 1;
% a = load('results/JMR/1-1Straight/polynome1-1Straight.mat');
% p_ini_v = a.p_ini_v;
% 
if move_dirichlet==1
    p_ini_v = [ones(l_spli*l_spli,1); ...
               ones(l_spli-2,1); ones(l_spli-2,1);1;1];
else
    p_ini_v = ones(l_spli*l_spli,1);
end

%%% Create shape parameters struct
param_sopt = struct;
param_sopt.friction = 0.8;
param_sopt.Fdes = Fdes;
param_sopt.Zdes = Zdes;
param_sopt.move_dirichlet = move_dirichlet;
fric = param_sopt.friction;

[sole.coor,dPFree_dp] = deformation_moveDiri(sole,p_ini_v,spl,move_dirichlet);


D = 3;
%sole.coor = deformation_moveDiri(sole,p_ini_v,spl,param_sopt.move_dirichlet); % Moving Dirichlet
% st1 = 1e-09;
% sole.coor(2,1) = sole.coor(2,1) + st1;
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
% for i=1:size(angleact_tot,2)
    angleact = angleact_tot(:,1);
    displ = displ_tot(:,1);
    contAngle = 1;
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
    FcOld3 = FcNew3;
    FcOld = reshape(FcOld3,3,sole.nFreeSurf);
    diplini = displ;
    angleactini = angleact;
    angleactOld = angleact;
    
    FtotZMPdes = [Fdes(:,1);Zdes(:,1);0];
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
    P0 = sole.coor';
    Pfree = P0(:,sole.nodesFreeSurf);
    Cs = sole.Cs;
    contAngle = 1;
    displ_first= [0;0;0];
    psi_first = 0;
    Ftot = FtotZMP(1:3);
    Z = FtotZMP(4:5);
    ind_cont3 = sort(([D*ind_cont-2; D*ind_cont-1; D*ind_cont]));
    ind_cont3 = reshape(ind_cont3,D*length(ind_cont),1); 
    Pfree_c = Pfree(:,ind_cont);
    Fc_mc = reshape(Fc_mc3,3,length(ind_cont));
    P_mc3 = PSurf3(ind_cont3);
    P_mc = reshape(P_mc3,3,length(ind_cont));
    theta = angleact(1);
    phi = angleact(2);
    psi = angleact(3);
    R = Rot(theta,phi,psi);
    Wc = zeros(3*length(ind_cont),3*length(ind_cont));
    for j = 1:length(ind_cont)
        for k = 1:length(ind_cont)
            Wc((3*j)-2:3*j,(3*k)-2:3*k) = R * Cs((3*ind_cont(j))-2:3*ind_cont(j),(3*ind_cont(k))-2:3*ind_cont(k)) * R';
        end
    end
    A = gradient_simu_A_contact(P_mc,Fc_mc,ind_cont,Wc,ind_slip,Pfree_c,param_sopt.friction,contAngle,displ_first,psi_first,angleact);
    der_C_s_dPfree = cell(3*sole.nFree,1);
    der_C_cc_dPfree = cell(3*sole.nFree,1);
    KNodesFreeSurf = reshape(sole.dof(sole.nodesFreeSurf,:)',[1,3*sole.nFreeSurf]);
    for j = 1:(3*sole.nTot)      
        der_C_s_dPfree{j} = sole.der_C{j}(KNodesFreeSurf,KNodesFreeSurf);
        der_C_cc_dPfree{j} = der_C_s_dPfree{j}(ind_cont3,ind_cont3);
    end  
    RBlock_c = [];
    for j=1:length(ind_cont)
        RBlock_c = blkdiag(RBlock_c,R);
    end
    psiini = 0;
    Rini = Rot(theta,phi,psiini);
    Rini_block_c = [];
    for j=1:length(ind_cont)
        Rini_block_c = blkdiag(Rini_block_c,Rini);
    end
    dF_dPfree = zeros(3*length(ind_cont),3*sole.nFree);
    %Theta_mat = zeros(3*length(ind_cont)+2*length(ind_slip),3);
    contact_surf = 1;
    contact = find(contact_surf==ind_cont,1);   
    if isempty(find(contact_surf==ind_cont,1))
        contact = 0;
    end
    Theta_mat = gradient_simu_Theta_contact(sole,ind_cont,der_C_cc_dPfree,ind_slip,Fc_mc3,R,Rini);
    inv_A = inv(A);
    inv_A_theta = inv_A * Theta_mat;
    [dOl_dPfree,dY_dPfree,dFc_dPfree,ddeltat_dPfree] = compute_gradient_step(sole,p_ini_v,displ,angleact,Zdes,fric,Pfree,PSurf3,PabsOld,Fc_mc3,Cs,dPFree_dp,ind_cont,ind_slip,displ_first,psi_first,contAngle);

    f = @(x) testderOl_Pfree_i1(sole,x);
    J = diff5points(f,[2,1],1e-12,coorini);
    (dOl_dPfree(:,4)-J)./dOl_dPfree(:,4)   