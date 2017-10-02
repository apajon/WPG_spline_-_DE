function der_delta1 = der_delta(deltat1,deltat2,Ft1,Ft2)
der_delta1 = zeros(2,2);
der_delta1(1,:)= [-(Ft1*deltat2*(deltat1 - deltat2))/(deltat1^2 + deltat2^2)^(3/2), (Ft1*deltat1*(deltat1 - deltat2))/(deltat1^2 + deltat2^2)^(3/2)];
der_delta1(2,:)= [-(Ft2*deltat2*(deltat1 - deltat2))/(deltat1^2 + deltat2^2)^(3/2), (Ft2*deltat1*(deltat1 - deltat2))/(deltat1^2 + deltat2^2)^(3/2)];