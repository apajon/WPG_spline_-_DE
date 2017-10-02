clear all
clc
% syms x1 x2 x3 y1 y2 y3 z1 z2 z3 z4 x4 y4 real
% syms mu lambda real
% M = [1 1 1 1; x1 x2 x3 x4; y1 y2 y3 y4; z1 z2 z3 z4];
% dPhi = inv(M)*[zeros(1,3);eye(3)];
% C(1:3,1:3) = lambda*ones(3,3)+2*mu*eye(3);
% C(4:6,4:6) = mu*eye(3);
% 
% Be = sym(zeros(6,12));
% Be([5,6,3],3:3:12) = dPhi';
% Be([1,4,5],1:3:10) = dPhi';
% Be([4,2,6],2:3:11) = dPhi';
% % stima = det(M)/6*(Be'*C*Be);
% %f = matlabFunction(stima,'file','stimaSymbOld.m');
% %stimaOld = stimaSymbOld(0.1,0.5,V);
% 
% 
% %R1 = sym(zeros(6,12));
% % R1 = Be * det(M);
% % stima = 1/(6*det(M)) * (R1'*C*R1);
% % der_st_dx1 = diff(stima,x1);
% % f = matlabFunction(der_st_dx1,'file','der_stima_dx1.m');
% 
% % A = simplify(det(M)*dPhi);
% % Be1 = sym(zeros(6,12));
% % Be1([1,4,5],1:3:10) = A';
% % Be1([4,2,6],2:3:11) = A';
% % Be1([5,6,3],3:3:12) = A';
% % %stima = 1/6*(1/det(M))*(Be1'*C*Be1);
% % der_stima_dx1_first = -1/6*(1/det(M))*trace(inv(M)*diff(M,x1))*(Be1'*C*Be1);
% % der_stima_dx1_second = 1/6*(1/det(M))*(diff(Be1,x1)'*C*Be1 + Be1'*C*diff(Be1,x1));
% % der_st_dx1 = der_stima_dx1_first + der_stima_dx1_second;
% % 
% % f = matlabFunction(der_st_dx1,'file','der_stima_dx1.m');
% 
% A = simplify(det(M)*dPhi);
% Be1 = sym(zeros(6,12));
% Be1([1,4,5],1:3:10) = A';
% Be1([4,2,6],2:3:11) = A';
% Be1([5,6,3],3:3:12) = A';
% %stima = 1/6*(1/det(M))*(Be1'*C*Be1);
% der_stima_dx2_first = -1/6*(1/det(M))*trace(inv(M)*diff(M,x2))*(Be1'*C*Be1);
% der_stima_dx2_second = 1/6*(1/det(M))*(diff(Be1,x2)'*C*Be1 + Be1'*C*diff(Be1,x2));
% der_st_dx2 = der_stima_dx2_first + der_stima_dx2_second;
% 
% f = matlabFunction(der_st_dx2,'file','der_stima_dx2.m');
% 
% A = simplify(det(M)*dPhi);
% Be1 = sym(zeros(6,12));
% Be1([1,4,5],1:3:10) = A';
% Be1([4,2,6],2:3:11) = A';
% Be1([5,6,3],3:3:12) = A';
% %stima = 1/6*(1/det(M))*(Be1'*C*Be1);
% der_stima_dx3_first = -1/6*(1/det(M))*trace(inv(M)*diff(M,x3))*(Be1'*C*Be1);
% der_stima_dx3_second = 1/6*(1/det(M))*(diff(Be1,x3)'*C*Be1 + Be1'*C*diff(Be1,x3));
% der_st_dx3 = der_stima_dx3_first + der_stima_dx3_second;
% 
% f = matlabFunction(der_st_dx3,'file','der_stima_dx3.m');
% 
% A = simplify(det(M)*dPhi);
% Be1 = sym(zeros(6,12));
% Be1([1,4,5],1:3:10) = A';
% Be1([4,2,6],2:3:11) = A';
% Be1([5,6,3],3:3:12) = A';
% %stima = 1/6*(1/det(M))*(Be1'*C*Be1);
% der_stima_dx4_first = -1/6*(1/det(M))*trace(inv(M)*diff(M,x4))*(Be1'*C*Be1);
% der_stima_dx4_second = 1/6*(1/det(M))*(diff(Be1,x4)'*C*Be1 + Be1'*C*diff(Be1,x4));
% der_st_dx4 = der_stima_dx4_first + der_stima_dx4_second;
% 
% f = matlabFunction(der_st_dx4,'file','der_stima_dx4.m');
% 
% A = simplify(det(M)*dPhi);
% Be1 = sym(zeros(6,12));
% Be1([1,4,5],1:3:10) = A';
% Be1([4,2,6],2:3:11) = A';
% Be1([5,6,3],3:3:12) = A';
% %stima = 1/6*(1/det(M))*(Be1'*C*Be1);
% der_stima_dy1_first = -1/6*(1/det(M))*trace(inv(M)*diff(M,y1))*(Be1'*C*Be1);
% der_stima_dy1_second = 1/6*(1/det(M))*(diff(Be1,y1)'*C*Be1 + Be1'*C*diff(Be1,y1));
% der_st_dy1 = der_stima_dy1_first + der_stima_dy1_second;
% 
% f = matlabFunction(der_st_dy1,'file','der_stima_dy1.m');
% 
% A = simplify(det(M)*dPhi);
% Be1 = sym(zeros(6,12));
% Be1([1,4,5],1:3:10) = A';
% Be1([4,2,6],2:3:11) = A';
% Be1([5,6,3],3:3:12) = A';
% %stima = 1/6*(1/det(M))*(Be1'*C*Be1);
% der_stima_dy2_first = -1/6*(1/det(M))*trace(inv(M)*diff(M,y2))*(Be1'*C*Be1);
% der_stima_dy2_second = 1/6*(1/det(M))*(diff(Be1,y2)'*C*Be1 + Be1'*C*diff(Be1,y2));
% der_st_dy2 = der_stima_dy2_first + der_stima_dy2_second;
% 
% f = matlabFunction(der_st_dy2,'file','der_stima_dy2.m');
% 
% A = simplify(det(M)*dPhi);
% Be1 = sym(zeros(6,12));
% Be1([1,4,5],1:3:10) = A';
% Be1([4,2,6],2:3:11) = A';
% Be1([5,6,3],3:3:12) = A';
% %stima = 1/6*(1/det(M))*(Be1'*C*Be1);
% der_stima_dy3_first = -1/6*(1/det(M))*trace(inv(M)*diff(M,y3))*(Be1'*C*Be1);
% der_stima_dy3_second = 1/6*(1/det(M))*(diff(Be1,y3)'*C*Be1 + Be1'*C*diff(Be1,y3));
% der_st_dy3 = der_stima_dy3_first + der_stima_dy3_second;
% 
% f = matlabFunction(der_st_dy3,'file','der_stima_dy3.m');
% 
% A = simplify(det(M)*dPhi);
% Be1 = sym(zeros(6,12));
% Be1([1,4,5],1:3:10) = A';
% Be1([4,2,6],2:3:11) = A';
% Be1([5,6,3],3:3:12) = A';
% %stima = 1/6*(1/det(M))*(Be1'*C*Be1);
% der_stima_dy4_first = -1/6*(1/det(M))*trace(inv(M)*diff(M,y4))*(Be1'*C*Be1);
% der_stima_dy4_second = 1/6*(1/det(M))*(diff(Be1,y4)'*C*Be1 + Be1'*C*diff(Be1,y4));
% der_st_dy4 = der_stima_dy4_first + der_stima_dy4_second;
% 
% f = matlabFunction(der_st_dy4,'file','der_stima_dy4.m');
% 
% A = simplify(det(M)*dPhi);
% Be1 = sym(zeros(6,12));
% Be1([1,4,5],1:3:10) = A';
% Be1([4,2,6],2:3:11) = A';
% Be1([5,6,3],3:3:12) = A';
% %stima = 1/6*(1/det(M))*(Be1'*C*Be1);
% der_stima_dz1_first = -1/6*(1/det(M))*trace(inv(M)*diff(M,z1))*(Be1'*C*Be1);
% der_stima_dz1_second = 1/6*(1/det(M))*(diff(Be1,z1)'*C*Be1 + Be1'*C*diff(Be1,z1));
% der_st_dz1 = der_stima_dz1_first + der_stima_dz1_second;
% 
% f = matlabFunction(der_st_dz1,'file','der_stima_dz1.m');
% 
% A = simplify(det(M)*dPhi);
% Be1 = sym(zeros(6,12));
% Be1([1,4,5],1:3:10) = A';
% Be1([4,2,6],2:3:11) = A';
% Be1([5,6,3],3:3:12) = A';
% %stima = 1/6*(1/det(M))*(Be1'*C*Be1);
% der_stima_dz2_first = -1/6*(1/det(M))*trace(inv(M)*diff(M,z2))*(Be1'*C*Be1);
% der_stima_dz2_second = 1/6*(1/det(M))*(diff(Be1,z2)'*C*Be1 + Be1'*C*diff(Be1,z2));
% der_st_dz2 = der_stima_dz2_first + der_stima_dz2_second;
% 
% f = matlabFunction(der_st_dz2,'file','der_stima_dz2.m');
% 
% A = simplify(det(M)*dPhi);
% Be1 = sym(zeros(6,12));
% Be1([1,4,5],1:3:10) = A';
% Be1([4,2,6],2:3:11) = A';
% Be1([5,6,3],3:3:12) = A';
% %stima = 1/6*(1/det(M))*(Be1'*C*Be1);
% der_stima_dz3_first = -1/6*(1/det(M))*trace(inv(M)*diff(M,z3))*(Be1'*C*Be1);
% der_stima_dz3_second = 1/6*(1/det(M))*(diff(Be1,z3)'*C*Be1 + Be1'*C*diff(Be1,z3));
% der_st_dz3 = der_stima_dz3_first + der_stima_dz3_second;
% 
% f = matlabFunction(der_st_dz3,'file','der_stima_dz3.m');
% 
% A = simplify(det(M)*dPhi);
% Be1 = sym(zeros(6,12));
% Be1([1,4,5],1:3:10) = A';
% Be1([4,2,6],2:3:11) = A';
% Be1([5,6,3],3:3:12) = A';
% %stima = 1/6*(1/det(M))*(Be1'*C*Be1);
% der_stima_dz4_first = -1/6*(1/det(M))*trace(inv(M)*diff(M,z4))*(Be1'*C*Be1);
% der_stima_dz4_second = 1/6*(1/det(M))*(diff(Be1,z4)'*C*Be1 + Be1'*C*diff(Be1,z4));
% der_st_dz4 = der_stima_dz4_first + der_stima_dz4_second;
% 
% f = matlabFunction(der_st_dz4,'file','der_stima_dz4.m');

% A = simplify(det(M)*dPhi);
% Be1 = sym(zeros(6,12));
% Be1([1,4,5],1:3:10) = A';
% Be1([4,2,6],2:3:11) = A';
% Be1([5,6,3],3:3:12) = A';
% stima = 1/6*(1/det(M))*(Be1'*C*Be1);
% der_stima_dx4_first = -1/6*(1/det(M))*trace(inv(M)*diff(M,x4))*(Be1'*C*Be1);
% der_stima_dx4_second = 1/6*(1/det(M))*(diff(Be1,x4)'*C*Be1 + Be1'*C*diff(Be1,x4));
% der_st_dx4 = der_stima_dx4_first + der_stima_dx4_second;
% % 
% % % diary off
% % % numop = feval(symengine, 'nops', der_t); 
% % % ops = cell(1,numop+1); 
% % % for K = 0 : numop 
% % % ops{K+1} = feval(symengine, 'op', der_t, K); 
% % % end
% % % der_p = diff(derCou,p);
% f = matlabFunction(der_st_dx4,'file','der_stima_dx4.m');
% 
% A = simplify(det(M)*dPhi);
% Be1 = sym(zeros(6,12));
% Be1([1,4,5],1:3:10) = A';
% Be1([4,2,6],2:3:11) = A';
% Be1([5,6,3],3:3:12) = A';
% stima = 1/6*(1/det(M))*(Be1'*C*Be1);
% der_stima_dy1_first = -1/6*(1/det(M))*trace(inv(M)*diff(M,y1))*(Be1'*C*Be1);
% der_stima_dy1_second = 1/6*(1/det(M))*(diff(Be1,y1)'*C*Be1 + Be1'*C*diff(Be1,y1));
% der_st_dy1 = der_stima_dy1_first + der_stima_dy1_second;
% % 
% % % diary off
% % % numop = feval(symengine, 'nops', der_t); 
% % % ops = cell(1,numop+1); 
% % % for K = 0 : numop 
% % % ops{K+1} = feval(symengine, 'op', der_t, K); 
% % % end
% % % der_p = diff(derCou,p);
% f = matlabFunction(der_st_dy1,'file','der_stima_dy1.m');
% 
% A = simplify(det(M)*dPhi);
% Be1 = sym(zeros(6,12));
% Be1([1,4,5],1:3:10) = A';
% Be1([4,2,6],2:3:11) = A';
% Be1([5,6,3],3:3:12) = A';
% stima = 1/6*(1/det(M))*(Be1'*C*Be1);
% der_stima_dy2_first = -1/6*(1/det(M))*trace(inv(M)*diff(M,y2))*(Be1'*C*Be1);
% der_stima_dy2_second = 1/6*(1/det(M))*(diff(Be1,y2)'*C*Be1 + Be1'*C*diff(Be1,y2));
% der_st_dy2 = der_stima_dy2_first + der_stima_dy2_second;
% % 
% % % diary off
% % % numop = feval(symengine, 'nops', der_t); 
% % % ops = cell(1,numop+1); 
% % % for K = 0 : numop 
% % % ops{K+1} = feval(symengine, 'op', der_t, K); 
% % % end
% % % der_p = diff(derCou,p);
% f = matlabFunction(der_st_dy2,'file','der_stima_dy2.m');
% 
% A = simplify(det(M)*dPhi);
% Be1 = sym(zeros(6,12));
% Be1([1,4,5],1:3:10) = A';
% Be1([4,2,6],2:3:11) = A';
% Be1([5,6,3],3:3:12) = A';
% stima = 1/6*(1/det(M))*(Be1'*C*Be1);
% der_stima_dy3_first = -1/6*(1/det(M))*trace(inv(M)*diff(M,y3))*(Be1'*C*Be1);
% der_stima_dy3_second = 1/6*(1/det(M))*(diff(Be1,y3)'*C*Be1 + Be1'*C*diff(Be1,y3));
% der_st_dy3 = der_stima_dy3_first + der_stima_dy3_second;
% % 
% % % diary off
% % % numop = feval(symengine, 'nops', der_t); 
% % % ops = cell(1,numop+1); 
% % % for K = 0 : numop 
% % % ops{K+1} = feval(symengine, 'op', der_t, K); 
% % % end
% % % der_p = diff(derCou,p);
% f = matlabFunction(der_st_dy3,'file','der_stima_dy3.m');
% 
% A = simplify(det(M)*dPhi);
% Be1 = sym(zeros(6,12));
% Be1([1,4,5],1:3:10) = A';
% Be1([4,2,6],2:3:11) = A';
% Be1([5,6,3],3:3:12) = A';
% stima = 1/6*(1/det(M))*(Be1'*C*Be1);
% der_stima_dy4_first = -1/6*(1/det(M))*trace(inv(M)*diff(M,y4))*(Be1'*C*Be1);
% der_stima_dy4_second = 1/6*(1/det(M))*(diff(Be1,y4)'*C*Be1 + Be1'*C*diff(Be1,y4));
% der_st_dy4 = der_stima_dy4_first + der_stima_dy4_second;
% % 
% % % diary off
% % % numop = feval(symengine, 'nops', der_t); 
% % % ops = cell(1,numop+1); 
% % % for K = 0 : numop 
% % % ops{K+1} = feval(symengine, 'op', der_t, K); 
% % % end
% % % der_p = diff(derCou,p);
% f = matlabFunction(der_st_dy4,'file','der_stima_dy4.m');
% 
% A = simplify(det(M)*dPhi);
% Be1 = sym(zeros(6,12));
% Be1([1,4,5],1:3:10) = A';
% Be1([4,2,6],2:3:11) = A';
% Be1([5,6,3],3:3:12) = A';
% stima = 1/6*(1/det(M))*(Be1'*C*Be1);
% der_stima_dz1_first = -1/6*(1/det(M))*trace(inv(M)*diff(M,z1))*(Be1'*C*Be1);
% der_stima_dz1_second = 1/6*(1/det(M))*(diff(Be1,z1)'*C*Be1 + Be1'*C*diff(Be1,z1));
% der_st_dz1 = der_stima_dz1_first + der_stima_dz1_second;
% % 
% % % diary off
% % % numop = feval(symengine, 'nops', der_t); 
% % % ops = cell(1,numop+1); 
% % % for K = 0 : numop 
% % % ops{K+1} = feval(symengine, 'op', der_t, K); 
% % % end
% % % der_p = diff(derCou,p);
% f = matlabFunction(der_st_dz1,'file','der_stima_dz1.m');
% 
% A = simplify(det(M)*dPhi);
% Be1 = sym(zeros(6,12));
% Be1([1,4,5],1:3:10) = A';
% Be1([4,2,6],2:3:11) = A';
% Be1([5,6,3],3:3:12) = A';
% stima = 1/6*(1/det(M))*(Be1'*C*Be1);
% der_stima_dz2_first = -1/6*(1/det(M))*trace(inv(M)*diff(M,z2))*(Be1'*C*Be1);
% der_stima_dz2_second = 1/6*(1/det(M))*(diff(Be1,z2)'*C*Be1 + Be1'*C*diff(Be1,z2));
% der_st_dz2 = der_stima_dz2_first + der_stima_dz2_second;
% % 
% % % diary off
% % % numop = feval(symengine, 'nops', der_t); 
% % % ops = cell(1,numop+1); 
% % % for K = 0 : numop 
% % % ops{K+1} = feval(symengine, 'op', der_t, K); 
% % % end
% % % der_p = diff(derCou,p);
% f = matlabFunction(der_st_dz2,'file','der_stima_dz2.m');
% 
% A = simplify(det(M)*dPhi);
% Be1 = sym(zeros(6,12));
% Be1([1,4,5],1:3:10) = A';
% Be1([4,2,6],2:3:11) = A';
% Be1([5,6,3],3:3:12) = A';
% stima = 1/6*(1/det(M))*(Be1'*C*Be1);
% der_stima_dz3_first = -1/6*(1/det(M))*trace(inv(M)*diff(M,z3))*(Be1'*C*Be1);
% der_stima_dz3_second = 1/6*(1/det(M))*(diff(Be1,z3)'*C*Be1 + Be1'*C*diff(Be1,z3));
% der_st_dz3 = der_stima_dz3_first + der_stima_dz3_second;
% % 
% % % diary off
% % % numop = feval(symengine, 'nops', der_t); 
% % % ops = cell(1,numop+1); 
% % % for K = 0 : numop 
% % % ops{K+1} = feval(symengine, 'op', der_t, K); 
% % % end
% % % der_p = diff(derCou,p);
% f = matlabFunction(der_st_dz3,'file','der_stima_dz3.m');
% 
% A = simplify(det(M)*dPhi);
% Be1 = sym(zeros(6,12));
% Be1([1,4,5],1:3:10) = A';
% Be1([4,2,6],2:3:11) = A';
% Be1([5,6,3],3:3:12) = A';
% stima = 1/6*(1/det(M))*(Be1'*C*Be1);
% der_stima_dz4_first = -1/6*(1/det(M))*trace(inv(M)*diff(M,z4))*(Be1'*C*Be1);
% der_stima_dz4_second = 1/6*(1/det(M))*(diff(Be1,z4)'*C*Be1 + Be1'*C*diff(Be1,z4));
% der_st_dz4 = der_stima_dz4_first + der_stima_dz4_second;
% % 
% % % diary off
% % % numop = feval(symengine, 'nops', der_t); 
% % % ops = cell(1,numop+1); 
% % % for K = 0 : numop 
% % % ops{K+1} = feval(symengine, 'op', der_t, K); 
% % % end
% % % der_p = diff(derCou,p);
% f = matlabFunction(der_st_dz4,'file','der_stima_dz4.m');

% stimaNew = stimaSymb(0.1,0.5,V);
V = [0.011321654348534  -0.038239115398949                   0;
     0.020339512021330  -0.044563429132578   0.030000000000000;
     0.016460559431850  -0.048549300362428                   0;
     0.030668309135475  -0.031393492702474   0.015308828067729];
% der_st_dx1 = der_stima_dx1(0.1,0.5,V);
% der_st_dx2 = der_stima_dx2(0.1,0.5,V);
% der_st_dx3 = der_stima_dx3(0.1,0.5,V);
% der_st_dx4 = der_stima_dx4(0.1,0.5,V);
% der_st_dy1 = der_stima_dy1(0.1,0.5,V);
% der_st_dy2 = der_stima_dy2(0.1,0.5,V);
% der_st_dy3 = der_stima_dy3(0.1,0.5,V);
% der_st_dy4 = der_stima_dy4(0.1,0.5,V);
% der_st_dz1 = der_stima_dz1(0.1,0.5,V);
% der_st_dz2 = der_stima_dz2(0.1,0.5,V);
% der_st_dz3 = der_stima_dz3(0.1,0.5,V);
% der_st_dz4 = der_stima_dz4(0.1,0.5,V);
lambda = 5.769230769230769e+05;
mu = 3.846153846153846e+05;
st1 = 1e-8;
V0 = V;
stima0 = stima(V0,lambda,mu);
%% x %%%
V(1,1) = V0(1,1) + st1;
stima1 = stima(V,lambda,mu);
G(:,:,1) = (stima1-stima0)/st1;
der_st_dx1 = der_stima_dx1(lambda,mu,V0);
der(:,:,1) = der_st_dx1;
% (G-der)./abs(der)
% %%%
V = V0;
V(2,1) = V0(2,1) + st1;
stima2 = stima(V,lambda,mu);
G(:,:,2) = (stima2-stima0)/st1;
der_st_dx2 = der_stima_dx2(lambda,mu,V0);
der(:,:,2) = der_st_dx2;
%%%
V = V0;
V(3,1) = V0(3,1) + st1;
stima3 = stima(V,lambda,mu);
G(:,:,3) = (stima3-stima0)/st1;
der_st_dx3 = der_stima_dx3(lambda,mu,V0);
der(:,:,3) = der_st_dx3;

%%%
V = V0;
V(4,1) = V0(4,1) + st1;
stima4 = stima(V,lambda,mu);
G(:,:,4) = (stima4-stima0)/st1;
der_st_dx4 = der_stima_dx4(lambda,mu,V0);
der(:,:,4) = der_st_dx4;


%% y %%%
%%%
V = V0;
V(1,2) = V0(1,2) + st1;
stima1 = stima(V,lambda,mu);
G(:,:,5) = (stima1-stima0)/st1;
der_st_dy1 = der_stima_dy1(lambda,mu,V0);
der(:,:,5) = der_st_dy1;
%%%
V = V0;
V(2,2) = V0(2,2) + st1;
stima2 = stima(V,lambda,mu);
G(:,:,6) = (stima2-stima0)/st1;
der_st_dy2 = der_stima_dy2(lambda,mu,V0);
der(:,:,6) = der_st_dy2;
%%%
V = V0;
V(3,2) = V0(3,2) + st1;
stima3 = stima(V,lambda,mu);
G(:,:,7) = (stima3-stima0)/st1;
der_st_dy3 = der_stima_dy3(lambda,mu,V0);
der(:,:,7) = der_st_dy3;
%%%
V = V0;
V(4,2) = V0(4,2) + st1;
stima4 = stima(V,lambda,mu);
G(:,:,8) = (stima4-stima0)/st1;
der_st_dy4 = der_stima_dy4(lambda,mu,V0);
der(:,:,8) = der_st_dy4;

%% z %%%
%%%
V = V0;
V(1,3) = V0(1,3) + st1;
stima1 = stima(V,lambda,mu);
G(:,:,9) = (stima1-stima0)/st1;
der_st_dz1 = der_stima_dz1(lambda,mu,V0);
der(:,:,9) = der_st_dz1;
%%%
V = V0;
V(2,3) = V0(2,3) + st1;
stima2 = stima(V,lambda,mu);
G(:,:,10) = (stima2-stima0)/st1;
der_st_dz2 = der_stima_dz2(lambda,mu,V0);
der(:,:,10) = der_st_dz2;
%%%
V = V0;
V(3,3) = V0(3,3) + st1;
stima3 = stima(V,lambda,mu);
G(:,:,11) = (stima3-stima0)/st1;
der_st_dz3 = der_stima_dz3(lambda,mu,V0);
der(:,:,11) = der_st_dz3;

%%%
V = V0;
V(4,3) = V0(4,3) + st1;
stima4 = stima(V,lambda,mu);
G(:,:,12) = (stima4-stima0)/st1;
der_st_dz4 = der_stima_dz4(lambda,mu,V0);
der(:,:,12) = der_st_dz4;

(G-der)./abs(der)
diff(Cou,qt) 