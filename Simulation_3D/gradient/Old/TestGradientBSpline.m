format long
clc
% clear all
close all

addpath ./input
addpath ./results/JMR/1-1Neoprene
addpath ./FEM
addpath ./iso2mesh
addpath ./geom3d/geom3d
addpath ./plane_line_intersect
addpath ./gradient
% addpath ./results/JMR/CartesianStiffness/Neopr

clear sole soleini
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %                       Sole FEM                               %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fname = 'semelle1.msh';
pname = 'input/semelle1 L=0.23, l=0.13, e=0.03 m new centre detailed/';
% %%% foot size %%%
l=0.13;
L=0.23;
e=0.03;
sole = soleFEM_newStiff(pname,fname,l,L,e);
coorini = sole.coor;
Young = 1000000;
Poisson = 0.3;
sole.setMaterial(Young,Poisson);
sole.stiffness();
% sole.derStiffSurface();
% K0 = sole.K;
% G1 = sole.der_K;
% %%%Gradient
% st1 = 1e-9;
% % G1 = sole.der_K;
% sole.coor(1,3) = coorini(1,3) + st1;
% sole.stiffness();
% K1 = sole.K;
% G = (K1-K0)/st1;
% E = (G - G1{3});
% E_rel = E/abs(G1{3});
% 
spline_res = -pi/2:pi/10:pi/2;
l_spli = length(spline_res);
spl = SplineClass(sole,spline_res);

move_dirichlet = 1;
% a = load('results/JMR/1-1Neoprene/polynome11Neoprene.mat');
% p_ini_v = a.p_ini_v;
% 
if move_dirichlet==1
    p_ini_v = [ones(l_spli*l_spli,1); ...
        ones(l_spli-2,1); ones(l_spli-2,1);1;1];
else
   p_ini_v = ones(l_spli*l_spli,1);
end
p0 = p_ini_v;
[coor0,G_spline_poly_0] = deformation_moveDiri(sole,p0,spl,move_dirichlet); % Moving Dirichlet
st = 1e-06;
p_ini_v(1)=p0(1)+st;
[coor1,G_spline_poly_1] = deformation_moveDiri(sole,p_ini_v,spl,move_dirichlet); % Moving Dirichlet
Gr = (reshape(coor1',3*sole.nTot,1)-reshape(coor0',3*sole.nTot,1))/st;
%Gr = (p_ini_v-p0)/st;
d = (G_spline_poly_0(:,1));
% (d-Gr)./(abs(Gr))
ffd = (coor1'-coor0')/st;
iffd = sum(ffd)~=0;
ga = reshape(G_spline_poly_0(:,1),3,570);
iga = sum(ga)~=0;
ffd(:,iffd)
ga(:,iga)