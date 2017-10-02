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
coorini = sole.coor;
Young = 1000000;
Poisson = 0.3;
sole.setMaterial(Young,Poisson);
nodeAnalised = 1;
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
G1 = full(sole.der_C{nodeAnalised});
iffd1 = G1~=0;
C0 = sole.C;
st1 = 1e-09;
sole.coor(sole.nodesFree(1),nodeAnalised) = coorini(sole.nodesFree(1),nodeAnalised) + st1;
sole.stiffness();
C1 = sole.C;
G = full(C1-C0)/st1;
iffd = G~=0;
%  [G1(1,:)' G(1,:)']
% % E = (G - G1{1});
% iffd = G~=0;
% E_rel = (G(iffd) - sole.der_C{D}(iffd))./abs(sole.der_C{D}(iffd));
