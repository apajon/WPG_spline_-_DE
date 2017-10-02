function [ZMP,F] = changeRef(ZMPx, ZMPy, Fx, Fy, Fz,backtoankle,fronttoankle,exttoankle,inttoankle,xpankle,ypankle,rightorleft,trajectory)
%%%%%%%%%% Foot %%%%%%%%%%
% ---------------------
% |                    |
% |                    |   
% |                    |  0.13
% |                    |    
% ---------------------
%          0.23
% taille du pied :
% distance en l’arrière du pied et la cheville : 0.1
% distance en l’avant du pied et la cheville : 0.13
% distance en l’extérieur du pied et la cheville : 0.075 (bord gauche dans notre cas)
% distance en l’intérieur du pied et la cheville : 0.055 (bord droit dans notre cas)
% Je t’ai mis en plus le tracé des courbes.
% backtoankle=0.098; %from back to ankle of foot
% fronttoankle=0.128; %from  front to ankle of foot
% exttoankle=0.076; %from exterior to ankle of foot
% inttoankle=0.054; %from interior to ankle of foot
% xpankle=1.257728354176543; %x coordinate of ankle position
% ypankle=-0.045000000000611; %y coordinate of ankle position
%rightorleft : 1 for right and -1 for left
%%%%%%%%%%%%%%%%%%%%%%%%%%
Lfooty = exttoankle+inttoankle;
Lfootx = backtoankle+fronttoankle;
traslx=xpankle-backtoankle;
trasly=ypankle-rightorleft*(Lfooty/2-inttoankle);
sidey = -(Lfooty/2):0.001:(Lfooty/2);
sidex1 = zeros(1,size(sidey,2));
sidex2 = zeros(1,size(sidey,2)) + Lfootx;
figure(4); clf;
subplot(1,2,1)
title('\fontsize{35}ZMP Trajectory')
axis([-0.05 0.25 -0.10 0.10])
set(gca, 'FontSize', 15)
xlabel('ZMP X [m]','fontsize',15);
ylabel('ZMP Y [m]','fontsize',15);
hold on
plot(sidex1,sidey,'LineWidth',2)
plot(sidex2,sidey,'LineWidth',2)
sidex = 0:0.001:Lfootx;
sidey1 = zeros(1,size(sidex,2)) - Lfooty/2;
sidey2 = zeros(1,size(sidex,2)) + Lfooty/2;
plot(sidex,sidey1,'LineWidth',2)
plot(sidex,sidey2,'LineWidth',2)
ZMPxNew = ZMPx-traslx;
ZMPyNew = ZMPy-trasly;
interval = 11:(length(ZMPxNew)-10);
if trajectory==2
    ZMPyNew = zeros(size(ZMPyNew));
end
plot(ZMPxNew(interval),ZMPyNew(interval),'LineWidth',4)
% plot(backtoankle,rightorleft*(Lfooty/2-inttoankle),'+')
hold off

ZMP = [ZMPxNew';ZMPyNew'];
if trajectory==2
    Fy = zeros(size(Fy));
end
F = [Fx';Fy';Fz'];

subplot(1,2,2)
set(gca, 'FontSize', 15)

h(1) = plot(F(1,interval),'k','LineWidth',4);
hold on;
h(2) = plot(F(2,interval),'r','LineWidth',4);
hold on;
h(3) = plot(F(3,interval),'b','LineWidth',4);
legend(h,'F_x','F_y','F_z','Location','north');
xlabel('Time Step','fontsize',15);
ylabel('ZMP Force [N]','fontsize',15);
axis([0 length(ZMPxNew(interval)) -100 1000])
title('\fontsize{35}ZMP Force')


% figure(2);
% % subplot(2,2,4)
% axis([-0.05 0.25 -0.10 0.10])
% set(gca, 'FontSize', 16)
% % title('\fontsize{55}ZMP Trajectory')
% title('\fontsize{16}ZMP Trajectory')
% xlabel('ZMP X [m]','fontsize',16);
% ylabel('ZMP Y [m]','fontsize',16);
% hold on
% plot(sidex1,sidey,'LineWidth',2)
% plot(sidex2,sidey,'LineWidth',2)
% plot(sidex,sidey1,'LineWidth',2)
% plot(sidex,sidey2,'LineWidth',2)
% 
% plot(ZMPxNew(interval),ZMPyNew(interval),'LineWidth',2)

figure(5);
subplot(2,2,3)
set(gca, 'FontSize', 16)
h(1) = plot(F(1,interval),'k','LineWidth',2);
hold on;
h(2) = plot(F(2,interval),'r','LineWidth',2);
hold on;
h(3) = plot(F(3,interval),'b','LineWidth',2);
legend(h,'F_x','F_y','F_z','Location','northeast');
xlabel('Time Step','fontsize',16);
ylabel('ZMP Force [N]','fontsize',16);
axis([0 length(ZMPxNew(interval)) -100 1400])
title('\fontsize{16}ZMP Force')