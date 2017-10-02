function N_ini = gradient_simu_N_ini_contact(ind_cont,Pfree_c3,Fc_mc3,angleact,psi_first,contAngle)
    t1 = [1 0 0];
    t2 = [0 1 0];
    N_ini = zeros(6,6);
    theta = angleact(1);
    phi = angleact(2);
    psi = angleact(3);    
    [R,Rtheta,Rphi,Rpsi] = Rot(theta,phi,psi);        
    Rpsiini = [cos(psi_first), -sin(psi_first), 0; sin(psi_first), cos(psi_first), 0; 0, 0, 1];            
    dRinidtheta = [0 0 0; 0 -sin(theta) -cos(theta); 0 cos(theta) -sin(theta)]*Rphi*Rpsiini;        
    dRinidphi = Rtheta*[-sin(phi) 0 cos(phi); 0 0 0; -cos(phi) 0 -sin(phi)]*Rpsiini;    
    PfreeFt1 = zeros(3,1);
    PfreeFt2 = zeros(3,1);
    PfreeFn = zeros(3,1);
    for i=1:length(ind_cont)
        PfreeFt1 = PfreeFt1 + Pfree_c3((3*i-2):3*i) * Fc_mc3(3*i-2);
        PfreeFt2 = PfreeFt2 + Pfree_c3((3*i-2):3*i) * Fc_mc3(3*i-1);
        PfreeFn = PfreeFn + Pfree_c3((3*i-2):3*i) * Fc_mc3(3*i);     
    end
    N_ini(4,4) = -t1 * dRinidtheta * PfreeFn;
    N_ini(4,5) = -t1 * dRinidphi * PfreeFn;
    N_ini(5,4) = t2 * dRinidtheta * PfreeFn;
    N_ini(5,5) = t2 * dRinidphi * PfreeFn;
    N_ini(6,4) = t1 * dRinidtheta * PfreeFt2 - t2 * dRinidtheta * PfreeFt1;
    N_ini(6,5) = t1 * dRinidphi * PfreeFt2 - t2 * dRinidphi * PfreeFt1;
end