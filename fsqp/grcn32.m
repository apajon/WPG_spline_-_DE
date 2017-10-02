function gradgj=grcn32(j,x,dummy,cd)

switch j
      case 1,
	 gradgj=[3.e0*x(1)*x(1), -6.e0, -4.e0];
      case 2,
	 gradgj=-1*ones(1,3);
end
return;
