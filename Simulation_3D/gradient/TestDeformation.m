format long
clc
%clear all
close all
addpath ./input
addpath ./results
addpath ./FEM
addpath ./iso2mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Sole FEM                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fname = 'semelle1.msh';
pname = 'input/semelle1 L=0.23, l=0.13, e=0.03 m new centre detailed/';
%%% foot size %%%
l=0.13;
L=0.23;
e=0.03;
sole = soleFEM_newStiff(pname,fname,l,L,e);
coorini = sole.coor;

coor = zeros(sole.nTot,3);
coor(:,1)= sole.coor(:,1) - sole.trasl(1);
coor(:,2) = sole.coor(:,2) - sole.trasl(2);
coor(:,3) = sole.coor(:,3) - sole.trasl(3) - sole.zpoles;
azimuth = zeros(1,size(coor,1));
elevation = zeros(1,size(coor,1));
r = zeros(1,size(coor,1));
for i=1:size(sole.coor,1);
    [azimuth(1,i),elevation(1,i),r(1,i)] = cart2sph(-coor(i,3),coor(i,2),coor(i,1));
end
spline_res = -pi/2:pi/10:pi/2;
spl = splineBasis([azimuth;elevation;r], spline_res, spline_res);

% With fix Dirichlet node
p_mid = ones(length(spline_res),length(spline_res));
p_mid(10,2)=3
p1 = ones(1,size(p_mid,1)+2);
p1(1:2) = 1.;
p1(3:end) = 1.;
% p1(end-1:end) = 1.;
p2 = ones(1,size(p_mid,1)+2);
p2(1:end) = 1.;
p2(end-1:end) = 1.;
p3 = ones(size(p_mid,1),1)+(ones(size(p_mid,1),1)-p1(3));
p3(1:end) = 1;
p4 = ones(size(p_mid,1),1)+(ones(size(p_mid,1),1)-p1(end-2));
p = [p1; p3 p_mid p4; p2];

% With moving Dirichlet node
% p = p_mid;
% pole_values = [p(1,1) p(1,2) p(1,size(p,1)-1) p(1,size(p,1)); p(end,1) p(end,2) p(end,size(p,1)-1) p(end,size(p,1))];

coornew = zeros(size(sole.coor));
coor_tmp = zeros(1,3);
for i=1:size(sole.coor,1)
    s1 = sum(sum(p(spl(i).idx_u, spl(i).idx_v).*spl(i).basisValues));
    [coor_tmp(1,1),coor_tmp(1,2),coor_tmp(1,3)] = sph2cart(spl(i).u,spl(i).v,s1*spl(i).r);
    coornew(i,:) = [coor_tmp(1,3),coor_tmp(1,2),-coor_tmp(1,1)];    
end

%%% translate because spl is in coornew-sole.trasl
coornew(:,1)= coornew(:,1) + sole.trasl(1);
coornew(:,2) = coornew(:,2) + sole.trasl(2);
coornew(:,3) = coornew(:,3) + sole.trasl(3) + sole.zpoles;

stressVM0 = zeros(sole.nTot,1);
plotsole(2,sole.elements_surf,coornew,stressVM0,[],[],[],[],-37.5,30);




