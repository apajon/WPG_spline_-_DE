function Fc_mc3 = testderFc_Pfree_i(sole,Pfree)
friction = 0.8;
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
    FcNew3 = zeros(D*m,1);
    displ_first= [0;0;0];
    psi_first = 0;
else
    Pfree = P0(:,sole.nodesFreeSurf);
    PabsOld = Pg(sole.nodesFreeSurf,:)';
    FcNew = FcFreeSurf'; 
    FcNew3 = reshape(FcNew,D*m,1);
    displ_first= [0;0;0];
    psi_first = 0;
end
FcOld3 = FcNew3;
FcOld = reshape(FcOld3,3,sole.nFreeSurf);
diplini = displ;
angleactini = angleact;
angleactOld = angleact;

[displ,angleact,FtotZMP,Fc_mc3_out,ind_cont_out,ind_slip_out,Mzmp,PSurf3] = Gauss(contAngle,friction,m,displ,angleact,FcNew3,Pfree,PabsOld,Cs);

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
end