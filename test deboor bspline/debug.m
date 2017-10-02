% close all;figure;hold on;
% plot(obj.A_xzmp*[7;0;15;15;0;5;10;-10;2;0;]+obj.B_xcom);
% plot(obj.A_xcom*[7;0;15;15;0;5;10;-10;5;6;]+obj.B_xcom,'g');
% plot(wpg_param.discretization(6)+1:sum(wpg_param.discretization(6:7)),obj.A_xzmp1*[7;0;15;15;0;5;10;2;0;+100]+obj.B_xzmp1,'-*m');
% plot(wpg_param.discretization(6)+2:sum(wpg_param.discretization(6:7)),obj.A_xzmp2*[7;0;15;15;0;5;10;2;0;+100]+obj.B_xzmp2,'-*r')
% % 0 10 7 0 15 15 0 5 10
% %%
% psa_abcdDSP=[...
%   ...
%    0.177668871189123 ;...
%    0.101038485438903 ;...
%    0.253120610680723 ;...
%    0.292589338889690 ;...
%   -0.000001804465378 ;...
%    2.173477797104588 ;...
%                    0 ;...
%    0.124462595332667 ;...
%    0.249462595332667 ;...
%   -0.000001804464425 ;...
%   ...
%   ...
%    0.046152253812035 ;...
%   -0.268171991317153 ;...
%    0.271351782980484 ;...
%   -0.016313862409985 ;...
%   -0.000003530177643 ;...
%    1.567847111224605 ;...
%                    0 ;...
%   -0.059343951996159 ;...
%    0.059343951996159 ;...
%   -0.000003530175111];
% 
% [0.125035382684594;...
%    0.043634062110479;...
%   -0.053455998280487;...
%    0.156795310665870;...
%    0.076856662424188;...
%   -0.431397732690450;...
%    0.217736633822889;...
%    0.083796049844527;...
%    0.426110477753154;...
%    0.250035382684594;...
%    0.043634062110479;...
%   -0.053455998280487;...
%    0.124921168158283;...
%    0.249921168158283;...
%    0.171129821613148;...
%    0.051918499805270;...
%    0.039369220636390;...
%   -0.073574677476811;...
%   -0.001780857934877;...
%   -0.378229097102691;...
%   -0.047918143024391;...
%    0.020574172207975;...
%   -0.769860159429530;...
%    0.047830057352015;...
%    0.020432656989128;...
%    0.918470172293313;...
%    0.073574677476811;...
%    0.001780857934877;...
%    0.378229097102691;...
%   -0.053507036729512;...
%    0.053507036729512;...
%   -0.049428622581360;...
%   -0.012658571848750;...
%    0.109915548635483]
% 
% %%
% figure(1);
% hold on;
% plot(zmp.A_xfcom*psa_abcdDSP(1:end/2-6)+zmp.B_xfcom)
% plot(zmp.A_yfcom*psa_abcdDSP(end/2+1:end-6)+zmp.B_yfcom,'g')
% %%
% figure(1);
% clf
% hold on;
% % plot(zmp.A_xzmp*psa_abcdDSP(1:end/2-6)+zmp.B_xzmp,zmp.A_yzmp*psa_abcdDSP(end/2+1:end-6)+zmp.B_yzmp)
% % plot(zmp.A_xazmp*psa_abcdDSP(1:end/2-6)+zmp.B_xazmp)
% % plot(zmp.A_yazmp*psa_abcdDSP(end/2+1:end-6)+zmp.B_yazmp,'g')
% xzmp=zmp.A_xzmp*psa_abcdDSP(1:end/2-6)+zmp.B_xzmp;
% yzmp=zmp.A_yzmp*psa_abcdDSP(end/2+1:end-6)+zmp.B_yzmp;
% xazmp=(xzmp(1:end-2)+xzmp(3:end)-2*xzmp(2:end-1))/(1/50)^2;
% yazmp=(yzmp(1:end-2)+yzmp(3:end)-2*yzmp(2:end-1))/(1/50)^2;
% plot(zmp.A_xazmp*psa_abcdDSP(1:end/2-6)+zmp.B_xazmp)
% plot(zmp.A_yazmp*psa_abcdDSP(end/2+1:end-6)+zmp.B_yazmp,'g')
% plot(xazmp,'--')
% plot(yazmp,'--g')
% %%
% figure(1);
% clf
% hold on;
% % plot(zmp.A_xzmp*psa_abcdDSP(1:end/2-6)+zmp.B_xzmp,zmp.A_yzmp*psa_abcdDSP(end/2+1:end-6)+zmp.B_yzmp)
% % plot(zmp.A_xazmp*psa_abcdDSP(1:end/2-6)+zmp.B_xazmp)
% % plot(zmp.A_yazmp*psa_abcdDSP(end/2+1:end-6)+zmp.B_yazmp,'g')
% xzmp=zmp.A_xcom*psa_abcdDSP(1:end/2-6)+zmp.B_xcom;
% yzmp=zmp.A_ycom*psa_abcdDSP(end/2+1:end-6)+zmp.B_ycom;
% xazmp=(xzmp(1:end-2)+xzmp(3:end)-2*xzmp(2:end-1))/(1/50)^2;
% yazmp=(yzmp(1:end-2)+yzmp(3:end)-2*yzmp(2:end-1))/(1/50)^2;
% plot((zmp.A_xfcom*psa_abcdDSP(1:end/2-6)+zmp.B_xfcom)/wpg_param.m)
% plot((zmp.A_yfcom*psa_abcdDSP(end/2+1:end-6)+zmp.B_yfcom)/wpg_param.m,'g')
% plot(xazmp,'--')
% plot(yazmp,'--g')
% %%
% xA_dif_zmp=obj.A_xzmp;
% % xA_dif_zmp=xA_dif_zmp(2:end);
% xB_dif_zmp=obj.B_xzmp;
% xB_dif_zmp=xB_dif_zmp(1:end)-trajectories_zmp.xpzmp;
% yA_dif_zmp=obj.A_yzmp;
% xH=xA_dif_zmp'*xA_dif_zmp;
% xB=xB_dif_zmp'*xA_dif_zmp;
% % xA_dif_zmp=xA_dif_zmp(2:end);
% yB_dif_zmp=obj.B_yzmp;
% yB_dif_zmp=yB_dif_zmp(1:end)-trajectories_zmp.ypzmp;
% yH=yA_dif_zmp'*yA_dif_zmp;
% yB=yB_dif_zmp'*yA_dif_zmp;
% H=[xH zeros(size(yH));zeros(size(xH)) yH];
% f=[xB yB]';
% opt = optimset('Algorithm', 'interior-point-convex','Display','iter','MaxIter',10000);
% psa_abcdDSP = quadprog(H,f,[],[],[],[],[],[],[],opt);
% 
% figure;
% clf
% axis equal
% hold on
% plot(trajectories_zmp.xpzmp,trajectories_zmp.ypzmp,'--')
% plot(obj.A_xzmp*psa_abcdDSP(1:end/2)+obj.B_xzmp,obj.A_yzmp*psa_abcdDSP(1+end/2:end)+obj.B_yzmp)
%%
xA_dif_zmp=zmp.A_xzmp;
xA_dif_zmp=xA_dif_zmp(2:end,:);
xB_dif_zmp=zmp.B_xzmp;
xB_dif_zmp=xB_dif_zmp(2:end)-trajectories_zmp.xpzmp;
xH=xA_dif_zmp'*xA_dif_zmp;
xB=xB_dif_zmp'*xA_dif_zmp;

yA_dif_zmp=zmp.A_yzmp;
yA_dif_zmp=yA_dif_zmp(2:end,:);
yB_dif_zmp=zmp.B_yzmp;
yB_dif_zmp=yB_dif_zmp(2:end)-trajectories_zmp.ypzmp;
yH=yA_dif_zmp'*yA_dif_zmp;
yB=yB_dif_zmp'*yA_dif_zmp;

H=[xH zeros(size(yH));zeros(size(xH)) yH];
f=[xB yB]';
opt = optimset('Algorithm', 'interior-point-convex','Display','iter','MaxIter',10000);
psa_abcdDSP = quadprog(H,f,[],[],[],[],[],[],[],opt);

figure;
clf
axis equal
title('ZMP projected on the ground')
xlabel('x(m)')
ylabel('y(m)')
hold on
plot(trajectories_zmp.xpzmp,trajectories_zmp.ypzmp,'--')
plot(zmp.A_xzmp*psa_abcdDSP(1:end/2)+zmp.B_xzmp,zmp.A_yzmp*psa_abcdDSP(1+end/2:end)+zmp.B_yzmp)
plot(zmp.A_xzmp*wpg_param.psa_abcd(1:length(wpg_param.psa_abcd)/2)+zmp.B_xzmp,zmp.A_yzmp*wpg_param.psa_abcd(length(wpg_param.psa_abcd)/2+1:end)+zmp.B_yzmp,'g')
hleg = legend('zmp with poly 5th','ZMP with least sqare diff','ZMP with bspline deboor','Location','EastOutside');
% Make the text of the legend italic and color it brown
set(hleg,'FontAngle','italic','TextColor',[.3,.2,.1])

figure;
clf
title('ZMP x axis')
ylabel('x(m)')
xlabel('t step')
hold on
plot(2:61,trajectories_zmp.xpzmp,'--')
plot(zmp.A_xzmp*psa_abcdDSP(1:end/2)+zmp.B_xzmp)
plot(zmp.A_xzmp*wpg_param.psa_abcd(1:length(wpg_param.psa_abcd)/2)+zmp.B_xzmp,'g')
hleg = legend('zmp with poly 5th','ZMP with least sqare diff','ZMP with bspline deboor','Location','EastOutside');
% Make the text of the legend italic and color it brown
set(hleg,'FontAngle','italic','TextColor',[.3,.2,.1])

figure;
clf
title('ZMP y axis')
ylabel('y(m)')
xlabel('t step')
hold on
plot(2:61,trajectories_zmp.ypzmp,'--')
plot(zmp.A_yzmp*psa_abcdDSP(1+end/2:end)+zmp.B_yzmp)
plot(zmp.A_yzmp*wpg_param.psa_abcd(length(wpg_param.psa_abcd)/2+1:end)+zmp.B_yzmp,'g')
hleg = legend('zmp with poly 5th','ZMP with least sqare diff','ZMP with bspline deboor','Location','EastOutside');
% Make the text of the legend italic and color it brown
set(hleg,'FontAngle','italic','TextColor',[.3,.2,.1])

%%
walking_param.psa_abcd(1:length(walking_param.psa_abcd)/2)'*zmp.A_xfcom2*walking_param.psa_abcd(1:length(walking_param.psa_abcd)/2)+zmp.B_xfcom2*walking_param.psa_abcd(1:length(walking_param.psa_abcd)/2)+zmp.C_xfcom2+...
    walking_param.psa_abcd(length(walking_param.psa_abcd)/2+1:end)'*zmp.A_xfcom2*walking_param.psa_abcd(length(walking_param.psa_abcd)/2+1:end)+zmp.B_xfcom2*walking_param.psa_abcd(length(walking_param.psa_abcd)/2+1:end)+zmp.C_xfcom2
%%
psa_abcdDSP(1:end/2)'*zmp.A_xfcom2*psa_abcdDSP(1:end/2)+zmp.B_xfcom2*psa_abcdDSP(1:end/2)+zmp.C_xfcom2+...
    psa_abcdDSP(end/2+1:end)'*zmp.A_xfcom2*psa_abcdDSP(end/2+1:end)+zmp.B_xfcom2*psa_abcdDSP(end/2+1:end)+zmp.C_xfcom2
%%
8.324233625157001e+03
2.142931506355848e+05
1.069100921239620e+02

%%
figure;
clf
title('ZMP acceleration y axis')
ylabel('acc y(m.s^-2)')
xlabel('t step')
hold on
plot((trajectories_zmp.ypzmp(3:end)+trajectories_zmp.ypzmp(1:end-2)-2*trajectories_zmp.ypzmp(2:end-1))/(1/50)^2)

% plot(zmp.A_xazmp*psa_abcdDSP(1:end/2)+zmp.B_xazmp,'--')
plot(zmp.A_yazmp*psa_abcdDSP(end/2+1:end)+zmp.B_yazmp,'--g')

% plot(zmp.A_xazmp*walking_param.psa_abcd(1:length(walking_param.psa_abcd)/2)+zmp.B_xazmp)
plot(zmp.A_yazmp*walking_param.psa_abcd(length(walking_param.psa_abcd)/2+1:end)+zmp.B_yazmp,'g')
hleg = legend('zmp with poly 5th','ZMP with least sqare diff','ZMP with bspline deboor','Location','EastOutside');
% Make the text of the legend italic and color it brown
set(hleg,'FontAngle','italic','TextColor',[.3,.2,.1])

%%
figure;
clf
title('COM x axis')
ylabel('y(m)')
xlabel('t step')
hold on
plot(2:61,trajectories_zmp.xpcom,'--')
% plot((trajectories_zmp.ypcom(3:end)+trajectories_zmp.ypcom(1:end-2)-2*trajectories_zmp.ypcom(2:end-1))/(1/50)^2)

% plot(zmp.A_xazmp*psa_abcdDSP(1:end/2)+zmp.B_xazmp,'--')
plot(zmp.A_xcom*psa_abcdDSP(1:end/2)+zmp.B_xcom,'--g')

% plot(zmp.A_xazmp*walking_param.psa_abcd(1:length(walking_param.psa_abcd)/2)+zmp.B_xazmp)
plot(zmp.A_xcom*walking_param.psa_abcd(1:length(walking_param.psa_abcd)/2)+zmp.B_xcom,'g')
hleg = legend('zmp with poly 5th','ZMP with least sqare diff','ZMP with bspline deboor','Location','EastOutside');
% Make the text of the legend italic and color it brown
set(hleg,'FontAngle','italic','TextColor',[.3,.2,.1])

figure;
clf
title('COM y axis')
ylabel('y(m)')
xlabel('t step')
hold on
plot(2:61,trajectories_zmp.ypcom,'--')
% plot((trajectories_zmp.ypcom(3:end)+trajectories_zmp.ypcom(1:end-2)-2*trajectories_zmp.ypcom(2:end-1))/(1/50)^2)

% plot(zmp.A_xazmp*psa_abcdDSP(1:end/2)+zmp.B_xazmp,'--')
plot(zmp.A_ycom*psa_abcdDSP(end/2+1:end)+zmp.B_ycom,'--g')

% plot(zmp.A_xazmp*walking_param.psa_abcd(1:length(walking_param.psa_abcd)/2)+zmp.B_xazmp)
plot(zmp.A_ycom*walking_param.psa_abcd(length(walking_param.psa_abcd)/2+1:end)+zmp.B_ycom,'g')
hleg = legend('zmp with poly 5th','ZMP with least sqare diff','ZMP with bspline deboor','Location','EastOutside');
% Make the text of the legend italic and color it brown
set(hleg,'FontAngle','italic','TextColor',[.3,.2,.1])

%%
figure;
clf
title('COM x axis')
ylabel('y(m)')
xlabel('t step')
hold on
xpcom=trajectories_zmp.xpcom;
plot((xpcom(3:end)+xpcom(1:end-2)-2*xpcom(2:end-1))/(1/50)^2,'--')

xpcom=zmp.A_xcom*psa_abcdDSP(1:end/2)+zmp.B_xcom;
plot((xpcom(3:end)+xpcom(1:end-2)-2*xpcom(2:end-1))/(1/50)^2,'--g')

xpcom=zmp.A_xcom*walking_param.psa_abcd(1:length(walking_param.psa_abcd)/2)+zmp.B_xcom;
plot((xpcom(3:end)+xpcom(1:end-2)-2*xpcom(2:end-1))/(1/50)^2,'g')
hleg = legend('zmp with poly 5th','ZMP with least sqare diff','ZMP with bspline deboor','Location','EastOutside');
% Make the text of the legend italic and color it brown
set(hleg,'FontAngle','italic','TextColor',[.3,.2,.1])

figure;
clf
title('COM y axis')
ylabel('y(m)')
xlabel('t step')
hold on
xpcom=trajectories_zmp.ypcom;
plot((xpcom(3:end)+xpcom(1:end-2)-2*xpcom(2:end-1))/(1/50)^2,'--')

xpcom=zmp.A_ycom*psa_abcdDSP(end/2+1:end)+zmp.B_ycom;
plot((xpcom(3:end)+xpcom(1:end-2)-2*xpcom(2:end-1))/(1/50)^2,'--g')

xpcom=zmp.A_ycom*walking_param.psa_abcd(length(walking_param.psa_abcd)/2+1:end)+zmp.B_ycom;
plot((xpcom(3:end)+xpcom(1:end-2)-2*xpcom(2:end-1))/(1/50)^2,'g')

hleg = legend('zmp with poly 5th','ZMP with least sqare diff','ZMP with bspline deboor','Location','EastOutside');
% Make the text of the legend italic and color it brown
set(hleg,'FontAngle','italic','TextColor',[.3,.2,.1])

figure;
clf
title('ZMP acceleration y axis')
ylabel('acc y(m.s^-2)')
xlabel('t step')
hold on
plot(trajectories_zmp.xpcom,trajectories_zmp.ypcom)

% plot(zmp.A_xazmp*psa_abcdDSP(1:end/2)+zmp.B_xazmp,'--')
plot(zmp.A_xcom*psa_abcdDSP(1:end/2)+zmp.B_xcom,zmp.A_ycom*psa_abcdDSP(end/2+1:end)+zmp.B_ycom,'--g')

% plot(zmp.A_xazmp*walking_param.psa_abcd(1:length(walking_param.psa_abcd)/2)+zmp.B_xazmp)
plot(zmp.A_xcom*walking_param.psa_abcd(1:length(walking_param.psa_abcd)/2)+zmp.B_xcom,zmp.A_ycom*walking_param.psa_abcd(length(walking_param.psa_abcd)/2+1:end)+zmp.B_ycom,'g')
hleg = legend('zmp with poly 5th','ZMP with least sqare diff','ZMP with bspline deboor','Location','EastOutside');
% Make the text of the legend italic and color it brown
set(hleg,'FontAngle','italic','TextColor',[.3,.2,.1])
%%
 switch (2)
    case 1:find(type_phase==0,1,'first')-1
        1
    case {find(type_phase==0,1,'last')+1:length(type_phase)}
        2
    otherwise
         0
 end
 
 %%
%  close all
figure;
clf
title('ZMP acceleration y axis')
ylabel('acc x(m.s^-2)')
xlabel('t step')
hold on
xzmp1=zmp.A_xzmp1*wpg_param.psa_abcdDSP(1:end/2)+zmp.B_xzmp1;
plot((xzmp1(3:end)+xzmp1(1:end-2)-2*xzmp1(2:end-1))/(1/200)^2,'b')
plot(zmp.A_xzmp1_acc*wpg_param.psa_abcdDSP(1:end/2)+zmp.B_xzmp1_acc,'g')
hleg = legend('zmp with poly 5th','ZMP with least sqare diff','ZMP with bspline deboor','Location','EastOutside');
% Make the text of the legend italic and color it brown
set(hleg,'FontAngle','italic','TextColor',[.3,.2,.1])

figure;
clf
title('ZMP acceleration y axis')
ylabel('acc y(m.s^-2)')
xlabel('t step')
hold on
yzmp1=zmp.A_yzmp1*wpg_param.psa_abcdDSP(end/2+1:end)+zmp.B_yzmp1;
plot((yzmp1(3:end)+yzmp1(1:end-2)-2*yzmp1(2:end-1))/(1/200)^2,'b')
plot(zmp.A_yzmp1_acc*wpg_param.psa_abcdDSP(end/2+1:end)+zmp.B_yzmp1_acc,'g')
hleg = legend('zmp with poly 5th','ZMP with least sqare diff','ZMP with bspline deboor','Location','EastOutside');
% Make the text of the legend italic and color it brown
set(hleg,'FontAngle','italic','TextColor',[.3,.2,.1])

%%
(wpg_param.psa_abcdDSP'*zmp.H*wpg_param.psa_abcdDSP+2*zmp.G*wpg_param.psa_abcdDSP+zmp.C)/wpg_param.nbpointdiscret
(wpg_param.psa_abcdDSP(1:end/2)'*wpg_param.epsilon*(zmp.A_xzmp1_acc2+zmp.A_xzmp2_acc2)*wpg_param.psa_abcdDSP(1:end/2)+2*wpg_param.epsilon*(zmp.B_xzmp1_acc2+zmp.B_xzmp2_acc2)*wpg_param.psa_abcdDSP(1:end/2)+wpg_param.epsilon*(zmp.C_xzmp1_acc2+zmp.C_xzmp2_acc2))/wpg_param.nbpointdiscret
(wpg_param.psa_abcdDSP(1:end/2)'*0.5*wpg_param.mu*(zmp.A_yt2_zmp1+zmp.A_yt2_zmp2)*wpg_param.psa_abcdDSP(1:end/2)+2*0.5*wpg_param.mu*(zmp.B_yt2_zmp1+zmp.B_yt2_zmp2)*wpg_param.psa_abcdDSP(1:end/2)+0.5*wpg_param.mu*(zmp.C_yt2_zmp1+zmp.C_yt2_zmp2))/wpg_param.nbpointdiscret