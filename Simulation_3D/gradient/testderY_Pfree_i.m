function angleactfin = testderY_Pfree_i(coorini,sole,Zdes,Fdes,fric,angleact,displ)
sole.coor = coorini;
sole.stiffness();
sole.stiffnessSurface();
contAngle = 1;

Pg = [];
FcFreeSurf = [];
D = 3;

m = sole.nFreeSurf;
Cs = sole.Cs;
P0 = sole.coor';
if contAngle==1
    Pfree = P0(:,sole.nodesFreeSurf);
    PabsOld(1,:) = P0(1,sole.nodesFreeSurf);
    PabsOld(2,:) = P0(2,sole.nodesFreeSurf);
    PabsOld(3,:) = P0(3,sole.nodesFreeSurf);    
    FSurfNew3 = zeros(D*m,1);
    displ_first= [0;0;0];
else
    Pfree = P0(:,sole.nodesFreeSurf);
    PabsOld = Pg(sole.nodesFreeSurf,:)';
    FcNew = FcFreeSurf'; 
    FSurfNew3 = reshape(FcNew,D*m,1);
    displ_first= [0;0;0];
end
FtotZMPdes = [Fdes;Zdes;0]; 
[displ,angleact,FtotZMP,Fc_mc3_out,PSurf3,ind_cont_out,ind_slip_out,Kcart,Bomega_out,J2,A_out,B_out] = GaussFtotZMP(contAngle,fric,m,displ,angleact,FtotZMPdes,FSurfNew3,Pfree,PabsOld,Cs,displ_first);
angleactfin = angleact;
end