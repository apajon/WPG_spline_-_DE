syms real x1 y1 z1 x2 y2 z2 x3 y3 z3 x4 y4 z4 lambda mu
V = sym('V',[4,3]);
V(1,1)=x1;
V(1,2)=y1;
V(1,3)=z1;
V(2,1)=x2;
V(2,2)=y2;
V(2,3)=z2;
V(3,1)=x3;
V(3,2)=y3;
V(3,3)=z3;
V(4,1)=x4;
V(4,2)=y4;
V(4,3)=z4;
M = [1 1 1 1; transpose(V)];
PhiGrad = M\[zeros(1,3);eye(3)];
Be = sym('Be',[6,12]);
for i=1:size(Be,1) 
    for j=1:size(Be,2) 
        Be(i,j)=0; 
    end
end
Be([1,4,5],1:3:10) = transpose(PhiGrad);
Be([4,2,6],2:3:11) = transpose(PhiGrad);
Be([5,6,3],3:3:12) = transpose(PhiGrad);
C(1:3,1:3) = lambda*ones(3,3)+2*mu*eye(3);
C(4:6,4:6) = mu*eye(3);
Be1 = det(M)*Be;

dBe1dx1 = diff(Be1,x1);