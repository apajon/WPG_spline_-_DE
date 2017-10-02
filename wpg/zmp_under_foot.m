function [time_ pzmp pcom fzmp gradpzmp gradfzmp] = zmp_under_foot(wpg_param,foot_step_wanted,nbpankle,time,trajectories_zmp,zpcom,zfzmp1,zmp,psa_abcd,discretization,discretization_,mg)
xpzmp_=[trajectories_zmp.xpzmp2;
    trajectories_zmp.xpzmp(end-(wpg_param.tss/2)*wpg_param.frequency+1:end);
    trajectories_zmp.xpzmp(1:sum(discretization(wpg_param.nbpolypi+wpg_param.nbpolyssp+wpg_param.nbpolydsp+wpg_param.nbpolyssp/2+1:wpg_param.nbpolypi+wpg_param.nbpolyssp+wpg_param.nbpolydsp+wpg_param.nbpolyssp)))+0.125;
    trajectories_zmp.xpzmp1+0.125];
ypzmp_=[trajectories_zmp.ypzmp2;
    trajectories_zmp.ypzmp(end-(wpg_param.tss/2)*wpg_param.frequency+1:end);
    -trajectories_zmp.ypzmp(1:sum(discretization(wpg_param.nbpolypi+wpg_param.nbpolyssp+wpg_param.nbpolydsp+wpg_param.nbpolyssp/2+1:wpg_param.nbpolypi+wpg_param.nbpolyssp+wpg_param.nbpolydsp+wpg_param.nbpolyssp)));
    -trajectories_zmp.ypzmp1];

xpcom_=[trajectories_zmp.xpcom(end-(wpg_param.tss/2)*wpg_param.frequency-wpg_param.tds*wpg_param.frequency+1+1:end);trajectories_zmp.xpcom(1:(wpg_param.tss/2)*wpg_param.frequency+wpg_param.tds*wpg_param.frequency)+0.125];
ypcom_=[trajectories_zmp.ypcom(end-(wpg_param.tss/2)*wpg_param.frequency-wpg_param.tds*wpg_param.frequency+1+1:end);-trajectories_zmp.ypcom(1:(wpg_param.tss/2)*wpg_param.frequency+wpg_param.tds*wpg_param.frequency)];
zpcom_=[zpcom(1)*ones(size(ypcom_))];

zfzmp2_=(1-zfzmp1(2:end))*mg;
zfzmp1_=(zfzmp1)*mg;
zfzmp=ones(wpg_param.tss*wpg_param.frequency,1)*mg;
zfzmp_=[zfzmp2_;zfzmp;zfzmp1_;];

xfzmp=-(zmp.A_xfcom*psa_abcd(1:length(psa_abcd)/2)+zmp.B_xfcom); % tu peux l'appeller Fx
yfzmp=-(zmp.A_yfcom*psa_abcd(length(psa_abcd)/2+1:end)+zmp.B_yfcom);
xfzmp_=[xfzmp(end-(wpg_param.tss/2)*wpg_param.frequency-wpg_param.tds*wpg_param.frequency+1+1:end);xfzmp(1:(wpg_param.tss/2)*wpg_param.frequency+wpg_param.tds*wpg_param.frequency)].*(zfzmp_/mg);
yfzmp_=[yfzmp(end-(wpg_param.tss/2)*wpg_param.frequency-wpg_param.tds*wpg_param.frequency+1+1:end);-yfzmp(1:(wpg_param.tss/2)*wpg_param.frequency+wpg_param.tds*wpg_param.frequency)].*(zfzmp_/mg);

pzmp=[xpzmp_ ypzmp_];
pcom=[xpcom_ ypcom_ zpcom_];
fzmp=[xfzmp_ yfzmp_ zfzmp_];

time_=(1:300)/200;


% load('trajectories.mat')
xApzmp_=[zmp.A_xzmp2;
    zmp.A_xzmp(end-(wpg_param.tss/2)*wpg_param.frequency+1:end,:) zeros(size(zmp.A_xzmp(end-(wpg_param.tss/2)*wpg_param.frequency+1:end,:),1),size(zmp.A_xzmp2,2)-length(psa_abcd)/2);
    zmp.A_xzmp(1:sum(discretization(wpg_param.nbpolypi+wpg_param.nbpolyssp+wpg_param.nbpolydsp+wpg_param.nbpolyssp/2+1:wpg_param.nbpolypi+wpg_param.nbpolyssp+wpg_param.nbpolydsp+wpg_param.nbpolyssp)),:) zeros(size(zmp.A_xzmp(1:sum(discretization(wpg_param.nbpolypi+wpg_param.nbpolyssp+wpg_param.nbpolydsp+wpg_param.nbpolyssp/2+1:wpg_param.nbpolypi+wpg_param.nbpolyssp+wpg_param.nbpolydsp+wpg_param.nbpolyssp)),:),1),size(zmp.A_xzmp2,2)-length(psa_abcd)/2);
    zmp.A_xzmp1];
yApzmp_=[zmp.A_yzmp2;
    zmp.A_yzmp(end-(wpg_param.tss/2)*wpg_param.frequency+1:end,:) zeros(size(zmp.A_yzmp(end-(wpg_param.tss/2)*wpg_param.frequency+1:end,:),1),size(zmp.A_yzmp2,2)-length(psa_abcd)/2);
    -zmp.A_yzmp(1:sum(discretization(wpg_param.nbpolypi+wpg_param.nbpolyssp+wpg_param.nbpolydsp+wpg_param.nbpolyssp/2+1:wpg_param.nbpolypi+wpg_param.nbpolyssp+wpg_param.nbpolydsp+wpg_param.nbpolyssp)),:) zeros(size(zmp.A_yzmp(1:sum(discretization(wpg_param.nbpolypi+wpg_param.nbpolyssp+wpg_param.nbpolydsp+wpg_param.nbpolyssp/2+1:wpg_param.nbpolypi+wpg_param.nbpolyssp+wpg_param.nbpolydsp+wpg_param.nbpolyssp)),:),1),size(zmp.A_yzmp2,2)-length(psa_abcd)/2);
    -zmp.A_yzmp1];
% zApzmp_=xApzmp_*0;

xgradpzmp=[xApzmp_ zeros(size(yApzmp_))];
ygradpzmp=[zeros(size(xApzmp_)) yApzmp_];
% zgradpzmp=xgradpzmp*0;
gradpzmp=[xgradpzmp;ygradpzmp];

xAfzmp_=diag(zfzmp_/mg)*[zmp.A_xfcom(end-(wpg_param.tss/2)*wpg_param.frequency-wpg_param.tds*wpg_param.frequency+1+1:end,:);zmp.A_xfcom(1:(wpg_param.tss/2)*wpg_param.frequency+wpg_param.tds*wpg_param.frequency,:)];
yAfzmp_=diag(zfzmp_/mg)*[zmp.A_yfcom(end-(wpg_param.tss/2)*wpg_param.frequency-wpg_param.tds*wpg_param.frequency+1+1:end,:);-zmp.A_yfcom(1:(wpg_param.tss/2)*wpg_param.frequency+wpg_param.tds*wpg_param.frequency,:)];
% zAfzmp_=xAfzmp_*0;

xgradfzmp=[-xAfzmp_ zeros(size(xAfzmp_,1),2*size(zmp.A_xzmp2,2)-length(psa_abcd)/2)];
ygradfzmp=[zeros(size(yAfzmp_,1),size(zmp.A_xzmp2,2)) -yAfzmp_ zeros(size(yAfzmp_,1),size(zmp.A_xzmp2,2)-length(psa_abcd)/2)];
% zgradfzmp=xgradfzmp*0;
gradfzmp=[xgradfzmp;ygradfzmp];

% %%
% trajectories=[time_' xpzmp_ ypzmp_ xpcom_ ypcom_ xfzmp yfzmp zfzmp];
% 
% zmpcom=fopen('exemple_trajectoire.txt','w');
% for i=1:size(trajectories,1)
%     fprintf(zmpcom,'%f %f %f %f %f %f %f %f\n',trajectories(i,:));
% end
% fclose(zmpcom);
% %%%%%%%%%%%
end