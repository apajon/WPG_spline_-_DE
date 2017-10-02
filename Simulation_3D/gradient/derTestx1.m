function D = derTestx1(lambda,mu,V)
M = [1 1 1 1; V'];
PhiGrad = M\[zeros(1,3);eye(3)];
Be = zeros(6,12);
Be([1,4,5],1:3:10) = PhiGrad';
Be([4,2,6],2:3:11) = PhiGrad';
Be([5,6,3],3:3:12) = PhiGrad';
C(1:3,1:3) = lambda*ones(3,3)+2*mu*eye(3);
C(4:6,4:6) = mu*eye(3);

y1=V(1,2);
z1=V(1,3);
y2=V(2,2);
z2=V(2,3);
y3=V(3,2);
z3=V(3,3);
y4=V(4,2);
z4=V(4,3);

dBe1dx1 = [ 0, 0, 0,       0,       0,       0,       0,       0,       0,       0,       0,       0;
 0, 0, 0,       0, z4 - z3,       0,       0, z2 - z4,       0,       0, z3 - z2,       0;
 0, 0, 0,       0,       0, y3 - y4,       0,       0, y4 - y2,       0,       0, y2 - y3;
 0, 0, 0, z4 - z3,       0,       0, z2 - z4,       0,       0, z3 - z2,       0,       0;
 0, 0, 0, y3 - y4,       0,       0, y4 - y2,       0,       0, y2 - y3,       0,       0;
 0, 0, 0,       0, y3 - y4, z4 - z3,       0, y4 - y2, z2 - z4,       0, y2 - y3, z3 - z2];

% dM = [ y3*z2 - y2*z3 + y2*z4 - y4*z2 - y3*z4 + y4*z3, 0, 0, 0;
%        y1*z3 - y3*z1 - y1*z4 + y4*z1 + y3*z4 - y4*z3, 0, 0, 0;
%        y2*z1 - y1*z2 + y1*z4 - y4*z1 - y2*z4 + y4*z2, 0, 0, 0;
%        y1*z2 - y2*z1 - y1*z3 + y3*z1 + y2*z3 - y3*z2, 0, 0, 0];
trace_DM = y3*z2 - y2*z3 + y2*z4 - y4*z2 - y3*z4 + y4*z3;

D = ((dBe1dx1'*C*Be + Be'*C*dBe1dx1)- trace_DM*(Be'*C*Be))/6;

end