function J = ffd(f,n,h,x0)
if length(n)==1
    r = 1:n;
else
    r = n;
end
J = {};
f0 = f(x0);
for i=r
    xp1 = x0;
    xp1(i) = xp1(i)+h;
    J{i} = (f(xp1) - f0)/h;
end
end