function gradgj=grcn32_bis(j,x,dummy,cd)

	 gradgj(1,:)=[3.e0*x(1)*x(1), -6.e0, -4.e0];
	 gradgj(2,:)=-1*ones(1,3);

return;
