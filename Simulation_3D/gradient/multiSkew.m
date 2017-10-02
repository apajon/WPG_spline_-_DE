function S = multiSkew(v)
n = length(v);
r1 = 1:3:n;
r2 = 2:3:n;
r3 = 3:3:n;
x = v(r1);
y = v(r2);
z = v(r3);
S = zeros(n,3);
S(r1,2) = -z;
S(r1,3) = y;
S(r2,1) = z;
S(r2,3) = -x;
S(r3,1) = -y;
S(r3,2) = x;
end