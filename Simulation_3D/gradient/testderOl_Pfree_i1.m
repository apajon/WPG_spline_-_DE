function displfin = testderOl_Pfree_i1(sole,Pfree)
Fdes = 1.0e+02 * [-0.052709310000000;0.030650450000000;1.279851840000000];
Zdes = [0.020625645823457;0.049298000000611];
fric = 0.8;
sole.coor = Pfree;
sole.stiffness();
sole.stiffnessSurface();
sole.derStiff();
% a = load('Pg.mat');
% Pg = a.Pg();
% b = load('FcFreeSurf.mat');
% FcFreeSurf = b.FcFreeSurf();
% a = load('displ.mat');
% displ = a.displ();
% b = load('angleact.mat');
% angleact = b.angleact();
% displ = 1.0e-03 * [0.115457331163844;0.006202241467888;-0.512942670322166];
% angleact = [0.000058485776918;-0.007589944959809;-0.000058477583626];
%ABe = prep_stressVonMises(sole);
%%% Test Gradient algo with given position and orientation
a = load('displ_tot.mat');
displ_tot = a.displ_tot(); 
b = load('angleact_tot.mat');
angleact_tot = b.angleact_tot();
Pg = [];
FcFreeSurf = [];
D = 3;

angleact = angleact_tot(:,1);
displ = displ_tot(:,1);
contAngle = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PsoleAngle.m
m = sole.nFreeSurf;
Cs = sole.Cs;
no_conv = [];
P0 = sole.coor';
if contAngle==1
    Pfree = P0(:,sole.nodesFreeSurf);
    PabsOld(1,:) = P0(1,sole.nodesFreeSurf);
    PabsOld(2,:) = P0(2,sole.nodesFreeSurf);
    PabsOld(3,:) = P0(3,sole.nodesFreeSurf);    
    FSurfNew3 = zeros(D*m,1);
    displ_first= [0;0;0];
    psi_first = 0;
else
    Pfree = P0(:,sole.nodesFreeSurf);
    PabsOld = Pg(sole.nodesFreeSurf,:)';
    FcNew = FcFreeSurf'; 
    FSurfNew3 = reshape(FcNew,D*m,1);
    displ_first= [0;0;0];
    psi_first = 0;
end

FtotZMPdes = [Fdes(:,1);Zdes(:,1);0];
[displ,angleact,FtotZMP,Fc_mc3_out,PSurf3,ind_cont_out,ind_slip_out,Kcart,Bomega_out,J2,A_out,B_out] = GaussFtotZMP(contAngle,fric,m,displ,angleact,FtotZMPdes,FSurfNew3,Pfree,PabsOld,Cs,displ_first);

if ind_cont_out(1)==0
    a = find(ind_cont_out==0,2);
    Fc_mc3 = Fc_mc3_out(1:D*(a(2)-1));
    ind_cont = (sort(ind_cont_out(1:(a(2)-1))+1))';      
else
    a = find(ind_cont_out==0,1);
    Fc_mc3 = Fc_mc3_out(1:D*(a-1));
    ind_cont = (sort(ind_cont_out(1:(a-1))+1))';              
end
if norm(ind_slip_out)>0
    if ind_slip_out(1)==0
        b = find(ind_slip_out==0,2);
        ind_slip = (sort(ind_slip_out(1:(b(2)-1))+1))';  
    else
        b = find(ind_slip_out==0,1);
        ind_slip = (sort(ind_slip_out(1:(b-1))+1))';      
    end
else
    ind_slip = [];
end
displfin = displ;

end