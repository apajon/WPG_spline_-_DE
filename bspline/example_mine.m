% Illustrates B-spline foot point calculation.
clear all 
close all
clc
% Copyright 2010 Levente Hunyadi

k = 4;
% knot sequence
t = [0 0 0 0 0.5 1 1 1 1];
f = 0:0.01:1; %  Is the discretization time
B = bspline_basismatrix(k,t,f);
% control points
D = [ 0.1993 0.4965 0.6671 0.7085 0.6809 ...
    ; 0.8377 0.8436 0.7617 0.6126 0.212];
% points on B-spline curve
%M = bspline_deboor(k,t,D,f);
C = (B*D')'; % interpolation by hand
% plot control points and spline
figure;
hold all;
plot(D(1,:), D(2,:), 'g');
%plot(M(1,:), M(2,:), 'b');
plot(C(1,:), C(2,:), 'r');
% plot(P(1,:), P(2,:), 'kx');
% plot(F(1,:), F(2,:), 'rx');
% for i = 1 : size(P,2)
%     line( ...
%         [P(1,i), F(1,i)], ...
%         [P(2,i), F(2,i)]);
% end
hold off;