function N_ini_surf = gradient_simu_N_ini_surf(ind_cont,Rini,P_mc3,Fc_mc3,angleact,psi_first,contAngle)

    N_ini_surf = zeros(6,6);
    theta = angleact(1);
    phi = angleact(2);
    psi = angleact(3);    
    [R,Rtheta,Rphi,Rpsi] = Rot(theta,phi,psi);        
    Rpsiini = [cos(psi_first), -sin(psi_first), 0; sin(psi_first), cos(psi_first), 0; 0, 0, 1];            
    dRinidtheta = [0 0 0; 0 -sin(theta) -cos(theta); 0 cos(theta) -sin(theta)]*Rphi*Rpsiini;        
    dRinidphi = Rtheta*[-sin(phi) 0 cos(phi); 0 0 0; -cos(phi) 0 -sin(phi)]*Rpsiini;    
    PfreeFt1 = 0;
    PfreeFt2 = 0;
    PfreeFn = 0;
    for i=1:length(ind_cont)
        PfreeFt1 = PfreeFt1 + P_mc3((3*i-2):3*i)*Fc_mc3((3*i)-2);
        PfreeFt2 = PfreeFt2 + P_mc3((3*i-2):3*i)*Fc_mc3((3*i)-1);
        PfreeFn = PfreeFn + P_mc3((3*i-2):3*i)*Fc_mc3((3*i));
    end
    N_ini_surf(4,4) = -dot(dRinidtheta(1,:),PfreeFn);
    N_ini_surf(4,5) = -dot(dRinidphi(1,:),PfreeFn);
    N_ini_surf(5,4) = dot(dRinidtheta(2,:),PfreeFn);
    N_ini_surf(5,5) = dot(dRinidphi(2,:),PfreeFn);    
    N_ini_surf(6,4) = dot(dRinidtheta(1,:),PfreeFt2) - dot(dRinidtheta(2,:),PfreeFt1);
    N_ini_surf(6,5) = dot(dRinidphi(1,:),PfreeFt2) - dot(dRinidphi(2,:),PfreeFt1);  