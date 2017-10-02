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

%sole.coor = deformation_moveDiri(sole,p_ini_v,spl,param_sopt.move_dirichlet); % Moving Dirichlet
st1 = 1e-09;
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
        FcNew3 = zeros(D*m,1);
        displ_first= [0;0;0];
        psi_first = 0;
    else
        Pfree = P0(:,sole.nodesFreeSurf);
        PabsOld = Pg(sole.nodesFreeSurf,:)';
        FcNew = FcFreeSurf'; 
        FcNew3 = reshape(FcNew,D*m,1);
        displ_first= [0;0;0];
        psi_first = 0;
    end
    FcOld3 = FcNew3;
    FcOld = reshape(FcOld3,3,sole.nFreeSurf);
    diplini = displ;
    angleactini = angleact;
    angleactOld = angleact;
    
    [displ,angleact,FtotZMP,Fc_mc3_out,ind_cont_out,ind_slip_out,Mzmp,PSurf3] = Gauss(contAngle,param_sopt.friction,m,displ,angleact,FcNew3,Pfree,PabsOld,Cs);

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
    PabsOld_mc = PabsOld(:,ind_cont);
    A = gradient_simu_A_contact(P_mc,Fc_mc,ind_cont,Wc,ind_slip,PabsOld_mc,param_sopt.friction,contAngle,displ_first,psi_first,angleact);
    der_C_s = cell(3*sole.nFreeSurf,1);
    KNodesFreeSurf = reshape(sole.dof(sole.nodesFreeSurf,:)',[1,3*sole.nFreeSurf]);
    for j = 1:(3*sole.nFreeSurf)      
        der_C_s{j} = sole.der_C{KNodesFreeSurf(j)}(KNodesFreeSurf,KNodesFreeSurf);
    end

    der_C_cc = cell(3*length(ind_cont),1);
    for j = 1:(3*length(ind_cont))
        der_C_cc{j} = der_C_s{ind_cont3(j)}(ind_cont3,ind_cont3);
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
    Theta_mat = gradient_simu_Theta_contact(ind_cont,der_C_cc,ind_slip,Fc_mc3,RBlock_c,Rini_block_c);
    G1 = inv(A) * Theta_mat;
    Gc = G1(1:3*length(ind_cont),1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ind_cont3 = sort([D*ind_cont-2 D*ind_cont-1 D*ind_cont]);
    ind_cont3 = reshape(ind_cont3,D*length(ind_cont),1); 
    theta = angleact(1);
    phi = angleact(2);
    psi = angleact(3);
    R = Rot(theta,phi,psi);
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
    
    Fc = FcFreeSurf(contact,:);
    Pc = Pg(sole.nodesFreeSurf(contact),:);
    
    stressVM = stressVonMises(sole,dP,ABe); % VonMises' Stress    
    plotsole(2,sole.elements_surf,Pg,stressVM,Pc,Fc,Z,Ftot,-127.5,30);

end
Fc_mc3_1 =[5.277001393952872;
  -5.770497547022657;
  16.330217660081061;
   0.934631080762176;
   1.763854768571514;
   5.256199423030246;
   5.848392716023073;
   0.683924627661593;
  13.881977605131153;
  -1.121459804048815;
  -7.433786549899944;
  17.782075171768877;
  -3.475016480096318;
  -1.090234739139082;
  10.709890125994381;
  -0.008187531125073;
   0.009147682606400;
   0.015345793549835;
  -0.923576732920039;
   1.745719616363874;
  17.761583750880863;
  -4.602214537726075;
   7.723057650907209;
  11.237910040523376;
  -6.528386166155969;
   5.052287512334161;
  10.318781044906686];

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
G = (Fc_mc3_1-Fc_mc3)/(1e-09);
[Gc G]
% G = reshape(G',size(G,1)*size(G,2),1);
displ1 = 1.0e-03 * [0.305168993283732;-0.326556692096206;-0.224672281864693];
displ0 = 1.0e-03 * [0.305168996514595;-0.326556691720616;-0.224672290409996];
G = (displ1-displ0)./1e-09;
% angleact1 = [-0.000298182558637;-0.009600218561009;-0.000326097207717];
E_rel = (G - G1)./abs(G1);
% angleact0 = [-0.000298471441569;-0.009599771727431;-0.000325495902487];