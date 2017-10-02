% function J = diff5points(f,coor,h,x0)
% % if length(n)==1
% %     r = 1:n;
% % else
% %     r = n;
% % end
% % J = {};
% % for i=rs
%     xm2 = x0;
%     xm2(coor(1),coor(2)) = xm2(coor(1),coor(2))-2*h;
%     xm1 = x0;
%     xm1(coor(1),coor(2)) = xm1(coor(1),coor(2))-h;
%     xp1 = x0;
%     xp1(coor(1),coor(2)) = xp1(coor(1),coor(2))+h;
%     xp2 = x0;
%     xp2(coor(1),coor(2)) = xp2(coor(1),coor(2))+2*h;
%     J = (-f(xp2) + 8*f(xp1)-8*f(xm1)+f(xm2))/(12*h);
% % end
% end
function J = diff5points(f,n,h,x0)
if length(n)==1
    r = 1:n;
else
    r = n;
end
J = {};
for i=r
    xm2 = x0;
    xm2(i) = xm2(i)-2*h;
    xm1 = x0;
    xm1(i) = xm1(i)-h;
    xp1 = x0;
    xp1(i) = xp1(i)+h;
    xp2 = x0;
    xp2(i) = xp2(i)+2*h;
    J{i} = (-f(xp2) + 8*f(xp1)-8*f(xm1)+f(xm2))/(12*h);
end
end