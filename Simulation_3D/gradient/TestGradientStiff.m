format long
clc
%clear all
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
coorini = sole.coor;
Young = 1000000;
Poisson = 0.3;
sole.setMaterial(Young,Poisson);
for pointAnalised = 1:1
for coordAnalised = 1:1
%%% Test Gradient Stiffness %%%
% sole.stiffness();
% sole.derStiff();
% K0 = sole.K;
% G1 = full(sole.der_K{nodeAnalised});
% iffd1 = G1~=0;
% st1 = 1e-6;
% sole.coor(sole.nodesFree(1),nodeAnalised) = coorini(sole.nodesFree(1),nodeAnalised) + st1;
% sole.stiffness();
% K1 = sole.K;
% G = full(K1-K0)/st1;
% iffd = G~=0;
% G_ana = (G1(iffd1));
% G_diff_ana = (G(iffd));
% if length(G_ana)<length(G_diff_ana)
%     G_ana = [G_ana;zeros(length(G_diff_ana)-length(G_ana),1)];
% end
% [G_ana G_diff_ana]
% G1(iffd)-G(iffd)
% E_rel = (G(iffd) - G1{nodeAnalised}(iffd))./abs(G1{nodeAnalised}(iffd));
%E_rel =(G - G1{nodeAnalised})./abs(G1{nodeAnalised});
%%% Test Gradient Compliance %%%
sole.stiffness();
sole.derStiff();
G1 = full(sole.der_A{3*(pointAnalised-1)+coordAnalised});
iffd1 = G1~=0;
A0 = sole.A;
st1 = 1e-09;
sole.coor(pointAnalised,coordAnalised) = coorini(pointAnalised,coordAnalised) + st1;
sole.stiffness();
A1 = sole.A;
G = full(A1-A0)/st1;
Dr = (G1-G)./G;
I0 = G==0;
Dr(I0) = 0;
[m1,mi1] = max(abs(Dr));
[m2,mi2] = max(m1);
disp(['coord ' num2str(coordAnalised)]);
disp(sprintf('relative max at %i,%i : %1.12d',mi1(mi2),mi2,Dr(mi1(mi2),mi2)));
disp(G(mi1(mi2),mi2));
disp(G1(mi1(mi2),mi2));
G0 = G1(I0);
disp(sprintf('abs max when G=0: %1.12d',max(abs(G0))));
end
end
%  [G1(1,:)' G(1,:)']
% % E = (G - G1{1});
% iffd = G~=0;
% E_rel = (G(iffd) - sole.der_C{D}(iffd))./abs(sole.der_C{D}(iffd));
