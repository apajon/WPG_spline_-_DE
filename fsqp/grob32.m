function gradfj=grob32(j,x,dummy,cd)

fa=2.e0*(x(1)+3.e0*x(2)+x(3));
fb=8.e0*(x(1)-x(2));
gradfj(1)=fa+fb;
gradfj(2)=fa*3.e0-fb;
gradfj(3)=fa;

return;
