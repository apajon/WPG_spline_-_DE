function [v6]=v6_n5_f(ti,i)
v6=[((-1).*(ti(i+4)+(-1).*ti(i)).^(-1).*(ti(i+5)+(-1).*ti(i)).^(-1).*ti(i).^5.*((-1).*ti(i)+ti(i+1)) ...
  .^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).*((-1).*ti(i)+ti(i+3)).^(-1)) (5.*(ti(i+4)+(-1).*ti(i)).^(-1) ...
  .*(ti(i+5)+(-1).*ti(i)).^(-1).*ti(i).^4.*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1) ...
  .*((-1).*ti(i)+ti(i+3)).^(-1)) ((-10).*(ti(i+4)+(-1).*ti(i)).^(-1).*(ti(i+5)+(-1).*ti(i)).^(-1) ...
  .*ti(i).^3.*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).*((-1).*ti(i)+ti(i+3)).^(-1)) ( ...
  10.*(ti(i+4)+(-1).*ti(i)).^(-1).*(ti(i+5)+(-1).*ti(i)).^(-1).*ti(i).^2.*((-1).*ti(i)+ti(i+1)).^( ...
  -1).*((-1).*ti(i)+ti(i+2)).^(-1).*((-1).*ti(i)+ti(i+3)).^(-1)) ((-5).*(ti(i+4)+(-1).*ti(i)).^(-1) ...
  .*(ti(i+5)+(-1).*ti(i)).^(-1).*ti(i).*((-1).*ti(i)+ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).*(( ...
  -1).*ti(i)+ti(i+3)).^(-1)) ((ti(i+4)+(-1).*ti(i)).^(-1).*(ti(i+5)+(-1).*ti(i)).^(-1).*((-1).*ti(i)+ ...
  ti(i+1)).^(-1).*((-1).*ti(i)+ti(i+2)).^(-1).*((-1).*ti(i)+ti(i+3)).^(-1))];
end