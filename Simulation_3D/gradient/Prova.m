Z1 = 1;
Z2 = 1;
Z3 = 2;
Z4 = 0;
Z5 = 3;
Z6 = 2;
Z7 = 5;
Z8 = 6;
Z9 = 7;
Z = [Z1 Z2 Z3; Z4 Z5 Z6; Z7 Z8 Z9];
Zt = Z';
[u,DE] = eig(Z);
u_max = u(:,1); % vector of matlab are normalized
[v,DV] = eig(Zt);
v_max = v(:,1);
%u_max = u_max/norm(u_max);
% dlambda = (v_max' *u_max)/(v_max' * u_max);
dlambda1 =  v_max'*[1 0 0; 0 0 0; 0 0 0] *u_max/(v_max' * u_max);
dlambda2 =  v_max'*[0 1 0; 0 0 0; 0 0 0] *u_max/(v_max' * u_max);
dlambda3 =  v_max'*[0 0 1; 0 0 0; 0 0 0] *u_max/(v_max' * u_max);
dlambda4 =  v_max'*[0 0 0; 1 0 0; 0 0 0] *u_max/(v_max' * u_max);
dlambda5 =  v_max'*[0 0 0; 0 1 0; 0 0 0] *u_max/(v_max' * u_max);
dlambda6 =  v_max'*[0 0 0; 0 0 1; 0 0 0] *u_max/(v_max' * u_max);
dlambda7 =  v_max'*[0 0 0; 0 0 0; 1 0 0] *u_max/(v_max' * u_max);
dlambda8 =  v_max'*[0 0 0; 0 0 0; 0 1 0] *u_max/(v_max' * u_max);
dlambda9 =  v_max'*[0 0 0; 0 0 0; 0 0 1] *u_max/(v_max' * u_max);

f = @(x) matrixEigenValue(x,Z2,Z3,Z4,Z5,Z6,Z7,Z8,Z9);
J1 = diff5points(f,1,1e-7,Z1);
(J1{1}-dlambda1)/J1{1}
f = @(x) matrixEigenValue(Z1,x,Z3,Z4,Z5,Z6,Z7,Z8,Z9);
J2 = diff5points(f,1,1e-7,Z2);
(J2{1}-dlambda2)/J2{1}
f = @(x) matrixEigenValue(Z1,Z2,x,Z4,Z5,Z6,Z7,Z8,Z9);
J3 = diff5points(f,1,1e-7,Z3);
(J3{1}-dlambda3)/J3{1}
f = @(x) matrixEigenValue(Z1,Z2,Z3,x,Z5,Z6,Z7,Z8,Z9);
J4 = diff5points(f,1,1e-7,Z4);
(J4{1}-dlambda4)/J4{1}
f = @(x) matrixEigenValue(Z1,Z2,Z3,Z4,x,Z6,Z7,Z8,Z9);
J5 = diff5points(f,1,1e-7,Z5);
(J5{1}-dlambda5)/J5{1}
f = @(x) matrixEigenValue(Z1,Z2,Z3,Z4,Z5,x,Z7,Z8,Z9);
J6 = diff5points(f,1,1e-7,Z6);
(J6{1}-dlambda6)/J6{1}
f = @(x) matrixEigenValue(Z1,Z2,Z3,Z4,Z5,Z6,x,Z8,Z9);
J7 = diff5points(f,1,1e-7,Z7);
(J7{1}-dlambda7)/J7{1}
f = @(x) matrixEigenValue(Z1,Z2,Z3,Z4,Z5,Z6,Z7,x,Z9);
J8 = diff5points(f,1,1e-7,Z8);
(J8{1}-dlambda8)/J8{1}
f = @(x) matrixEigenValue(Z1,Z2,Z3,Z4,Z5,Z6,Z7,Z8,x);
J9 = diff5points(f,1,1e-7,Z9);
(J9{1}-dlambda9)/J9{1}