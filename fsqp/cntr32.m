function gj=cntr32(j,x,cd)

switch j
      case 1
	 gj=x(1)^3-6.e0*x(2)-4.e0*x(3)+3.e0;
      case 2
	 gj=1.e0-x(1)-x(2)-x(3);
end

return;
