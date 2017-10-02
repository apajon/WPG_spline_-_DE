function RARt = blockMult(A,R)
     m3 = size(A,1); 
     m = m3/3;
%     I = reshape(1:(m3*m3),m3,m3);
%     J = reshape(I(:,[1:3:m3 2:3:m3 3:3:m3]),m*m3,3);
%     i = reshape(reshape(1:m3,m,3)',1,m3);
%     iJ = I(:,i);
%     tmp = reshape(R*reshape(A,3,m*m3),m3,m3);
%     tmp2 = tmp(J)*R';
%     RARt = tmp2(iJ);
    tmp = reshape(R*reshape(A,3,m*m3),m3,m3);
    tmp2 = tmp';
    RARt = reshape(R*reshape(tmp2,3,m*m3),m3,m3);
end