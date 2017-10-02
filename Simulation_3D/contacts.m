function simu = contacts(sole,simu,friction,contAngle)

if contAngle==1  
    simu.PSurfOld = simu.PfreeSurf;   
    simu.FSurf3 = zeros(3*sole.nFreeSurf,1);
else
    simu.PSurfOld = simu.P(sole.nodesFreeSurf,:)';
    simu.FSurf3 = simu.FSurf3; %%% no if not there are some problems for convergence 
    %simu.FSurf3 = zeros(3*sole.nFreeSurf,1);
    PSurfold = simu.P(sole.nodesFreeSurf,:)';
end
ind_contOld = simu.ind_cont;
ind_slipOld = simu.ind_slip;
Fold3 = simu.Fc3;
% simu.FSurf3 = zeros(3*sole.nFreeSurf,1);
diplini = simu.displ;
angleactini = simu.angleact;
Kcartold = simu.Kcart;
PSurfold3 = simu.PSurf3;
FtotZMPOld = simu.FtotZMP;
% a = load('PSurfOld.mat');
% simu.PSurfOld = a.PSurfOld;
% a = load('diplini.mat');
% simu.displ = a.diplini;
% a = load('angleactini.mat');
% simu.angleact = a.angleactini;

% save 'diplini.mat' diplini
% save 'angleactini.mat' angleactini
FtotZMPdes = [simu.Fdes;simu.Zdes;0]; 
[simu.displ,simu.angleact,simu.FtotZMP,Fc3_out,simu.PSurf3,ind_cont_out,ind_slip_out,simu.Kcart] = GaussFtotZMP(contAngle,friction,sole.nFreeSurf,simu.displ,simu.angleact,FtotZMPdes,simu.FSurf3,simu.PfreeSurf,simu.PSurfOld,sole.Cs,simu.displ_first);
% if norm(simu.Kcart)==0
%     simu.Kcart
% end
critZMP = norm(FtotZMPdes-simu.FtotZMP);
%norm(simu.Kcart)
if norm(simu.Kcart)==0
    simu.displ = diplini;
    simu.angleact = angleactini;
    simu.no_conv = contAngle;
    simu.Fc3 = Fold3;
    simu.ind_cont = ind_contOld;
    simu.ind_slip = ind_slipOld;
    simu.Kcart = Kcartold;
    simu.PSurf3 = PSurfold3;
    simu.FtotZMP = FtotZMPOld;    
else
    if norm(critZMP)>1e-03
        simu.displ = diplini;
        simu.angleact = angleactini;
        simu.no_conv = contAngle;
        simu.Fc3 = Fold3;
        simu.ind_cont = ind_contOld;
        simu.ind_slip = ind_slipOld;
        simu.Kcart = Kcartold;
        simu.PSurf3 = PSurfold3;
        simu.FtotZMP = FtotZMPOld;
    else
        if ind_cont_out(1)==0
            a = find(ind_cont_out==0,2);
            simu.Fc3 = Fc3_out(1:3*(a(2)-1));
            simu.ind_cont = (sort(ind_cont_out(1:(a(2)-1))+1))';      
        else
            a = find(ind_cont_out==0,1);
            simu.Fc3 = Fc3_out(1:3*(a-1));
            simu.ind_cont = (sort(ind_cont_out(1:(a-1))+1))';              
        end
        if norm(ind_slip_out)>0
            if ind_slip_out(1)==0
                b = find(ind_slip_out==0,2);
                simu.ind_slip = (sort(ind_slip_out(1:(b(2)-1))+1))';  
            else
                b = find(ind_slip_out==0,1);
                simu.ind_slip = (sort(ind_slip_out(1:(b-1))+1))';      
            end
        else
            simu.ind_slip = [];
        end
    end
end
simu.updateContact();
simu.updatePosition(sole);
simu.Ftot = simu.FtotZMP(1:3);
simu.Z = simu.FtotZMP(4:5);
% A_alg = A_out(1:(3*length(ind_cont)+2*length(ind_slip))*(3*length(ind_cont)+2*length(ind_slip)));
% A_alg = reshape(A_alg,3*length(ind_cont)+2*length(ind_slip),3*length(ind_cont)+2*length(ind_slip));
% B_alg = B_out(1:(3*simu.nc+2*simu.ns)*6);
% B_alg = reshape(B_alg,3*simu.nc+2*simu.ns,6);
