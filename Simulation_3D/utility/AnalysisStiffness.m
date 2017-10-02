% format long
% clc
% clear all
% close all
% 
% % addpath ./results/CMAME/CartesianStiffness/ButylRubber
% % addpath ./results/CMAME/CartesianStiffness/Neoprene
% % addpath ./results/CMAME/CartesianStiffness/NeopreneStraight
% % addpath ./results/CMAME/CartesianStiffness/1-1ButylRubber
% % addpath ./results/CMAME/CartesianStiffness/1-1Neoprene
% % addpath ./results/CMAME/CartesianStiffness/1-1NeoprenenoDir
% % addpath ./results/CMAME/CartesianStiffness/1-1NeopreneStraight
% % addpath ./results/CMAME/CartesianStiffness/1-1000Neoprene
% 
% a = load('./results/CMAME/CartesianStiffness/Neoprene/Kcart_vector.mat');
% Kcart_NE_vec = a.Kcart_vector;
% cost_rot_tot_NE = zeros(1,size(Kcart_NE_vec,1));
% cost_lin_tot_NE = zeros(1,size(Kcart_NE_vec,1));
% a = load('./results/CMAME/CartesianStiffness/ButylRubber/Kcart_vector.mat');
% Kcart_BR_vec = a.Kcart_vector;
% cost_rot_tot_BR = zeros(1,size(Kcart_NE_vec,1));
% cost_lin_tot_BR = zeros(1,size(Kcart_NE_vec,1));
% a = load('./results/CMAME/CartesianStiffness/Straight/Kcart_vector.mat');
% Kcart_Straight_vec = a.Kcart_vector;
% cost_rot_tot_Straight = zeros(1,size(Kcart_NE_vec,1));
% cost_lin_tot_Straight = zeros(1,size(Kcart_NE_vec,1));
% a = load('./results/CMAME/CartesianStiffness/1-1Neoprene/Kcart_vector.mat');
% Kcart_11NE_vec = a.Kcart_vector;
% cost_rot_tot_11NE = zeros(1,size(Kcart_NE_vec,1));
% cost_lin_tot_11NE = zeros(1,size(Kcart_NE_vec,1));
% a = load('./results/CMAME/CartesianStiffness/1-1ButylRubber/Kcart_vector.mat');
% Kcart_11BR_vec = a.Kcart_vector;
% cost_rot_tot_11BR = zeros(1,size(Kcart_NE_vec,1));
% cost_lin_tot_11BR = zeros(1,size(Kcart_NE_vec,1));
% a = load('./results/CMAME/CartesianStiffness/1-1Straight/Kcart_vector.mat');
% Kcart_11Straight_vec = a.Kcart_vector;
% cost_rot_tot_11Straight = zeros(1,size(Kcart_NE_vec,1));
% cost_lin_tot_11Straight = zeros(1,size(Kcart_NE_vec,1));
% a = load('./results/CMAME/CartesianStiffness/1-1000Neoprene/Kcart_vector.mat');
% Kcart_11000NE_vec = a.Kcart_vector;
% cost_rot_tot_11000NE = zeros(1,size(Kcart_NE_vec,1));
% cost_lin_tot_11000NE = zeros(1,size(Kcart_NE_vec,1));
% a = load('./results/CMAME/CartesianStiffness/1-1NoDir/Kcart_vector.mat');
% Kcart_noDir_vec = a.Kcart_vector;
% cost_rot_tot_noDir = zeros(1,size(Kcart_NE_vec,1));
% cost_lin_tot_noDir = zeros(1,size(Kcart_NE_vec,1));
% for i=1:size(Kcart_NE_vec,1)
%     Kcart_NE = Kcart_NE_vec{i};
%     cost_rot_tot_NE(i) = sqrt(min(eig(Kcart_NE(4:6,4:6)'*Kcart_NE(4:6,4:6))));
%     cost_lin_tot_NE(i) = norm(Kcart_NE(1:3,1:3),2);
%     Kcart_BR = Kcart_BR_vec{i};
%     cost_rot_tot_BR(i) = sqrt(min(eig(Kcart_BR(4:6,4:6)'*Kcart_BR(4:6,4:6))));
%     cost_lin_tot_BR(i) = norm(Kcart_BR(1:3,1:3),2);
%     Kcart_Straight = Kcart_Straight_vec{i};
%     cost_rot_tot_Straight(i) = sqrt(min(eig(Kcart_Straight(4:6,4:6)'*Kcart_Straight(4:6,4:6))));
%     cost_lin_tot_Straight(i) = norm(Kcart_Straight(1:3,1:3),2);
%     Kcart_11NE = Kcart_11NE_vec{i};
%     cost_rot_tot_11NE(i) = sqrt(min(eig(Kcart_11NE(4:6,4:6)'*Kcart_11NE(4:6,4:6))));
%     cost_lin_tot_11NE(i) = norm(Kcart_11NE(1:3,1:3),2);    
%     Kcart_11BR = Kcart_11BR_vec{i};
%     cost_rot_tot_11BR(i) = sqrt(min(eig(Kcart_11BR(4:6,4:6)'*Kcart_11BR(4:6,4:6))));
%     cost_lin_tot_11BR(i) = norm(Kcart_11BR(1:3,1:3),2);        
%     Kcart_11Straight = Kcart_11Straight_vec{i};
%     cost_rot_tot_11Straight(i) = sqrt(min(eig(Kcart_11Straight(4:6,4:6)'*Kcart_11Straight(4:6,4:6))));
%     cost_lin_tot_11Straight(i) = norm(Kcart_11Straight(1:3,1:3),2);
%     Kcart_11000NE = Kcart_11000NE_vec{i};
%     cost_rot_tot_11000NE(i) = sqrt(min(eig(Kcart_11000NE(4:6,4:6)'*Kcart_11000NE(4:6,4:6))));
%     cost_lin_tot_11000NE(i) = norm(Kcart_11000NE(1:3,1:3),2);
%     Kcart_noDir = Kcart_noDir_vec{i};
%     cost_rot_tot_noDir(i) = sqrt(min(eig(Kcart_noDir(4:6,4:6)'*Kcart_noDir(4:6,4:6))));
%     cost_lin_tot_noDir(i) = norm(Kcart_noDir(1:3,1:3),2);        
% end
% cost_lin_ini_NE = 0;
% cost_lin_ini_BR = 0;
% cost_lin_ini_Straight = 0;
% cost_lin_fin_11NE = 0;
% cost_lin_fin_11BR = 0;
% cost_lin_fin_11Straight = 0;
% cost_lin_fin_11000NE = 0;
% cost_lin_fin_noDir = 0;
% for i=1:(size(Kcart_NE_vec,1)/10)
%     cost_lin_ini_NE = cost_lin_ini_NE + cost_lin_tot_NE(i);
%     cost_lin_ini_BR = cost_lin_ini_BR + cost_lin_tot_BR(i);
%     cost_lin_ini_Straight = cost_lin_ini_Straight + cost_lin_tot_Straight(i);
%     cost_lin_fin_11NE = cost_lin_fin_11NE + cost_lin_tot_11NE(i);
%     cost_lin_fin_11BR = cost_lin_fin_11BR + cost_lin_tot_11BR(i);
%     cost_lin_fin_11Straight = cost_lin_fin_11Straight + cost_lin_tot_11Straight(i);
%     cost_lin_fin_11000NE = cost_lin_fin_11000NE + cost_lin_tot_11000NE(i);
%     cost_lin_fin_noDir = cost_lin_fin_noDir + cost_lin_tot_noDir(i);
% end
% cost_rot_ini_NE = 0;
% cost_rot_ini_BR = 0;
% cost_rot_ini_Straight = 0;
% cost_rot_fin_11NE = 0;
% cost_rot_fin_11BR = 0;
% cost_rot_fin_11Straight = 0;
% cost_rot_fin_11000NE = 0;
% cost_rot_fin_noDir = 0;
% for i=1:size(Kcart_NE_vec,1)
%     cost_rot_ini_NE = cost_rot_ini_NE + cost_rot_tot_NE(i);
%     cost_rot_ini_BR = cost_rot_ini_BR + cost_rot_tot_BR(i);
%     cost_rot_ini_Straight = cost_rot_ini_Straight + cost_rot_tot_Straight(i);
%     cost_rot_fin_11NE = cost_rot_fin_11NE + cost_rot_tot_11NE(i);
%     cost_rot_fin_11BR = cost_rot_fin_11BR + cost_rot_tot_11BR(i);
%     cost_rot_fin_11Straight = cost_rot_fin_11Straight + cost_rot_tot_11Straight(i);
%     cost_rot_fin_11000NE = cost_rot_fin_11000NE + cost_rot_tot_11000NE(i);
%     cost_rot_fin_noDir = cost_rot_fin_noDir + cost_rot_tot_noDir(i);
% end
% cost_ini_11NE = cost_lin_ini_NE - cost_rot_ini_NE;
% cost_ini_11000NE = cost_lin_ini_NE - 1000*cost_rot_ini_NE;
% cost_fin_11NE = cost_lin_fin_11NE - cost_rot_fin_11NE;
% cost_fin_11000NE = cost_lin_fin_11000NE - 1000*cost_rot_fin_11000NE;
% cost_red_perc_11NE = 100-(abs(cost_ini_11NE-cost_fin_11NE)/abs(cost_ini_11NE) * 100);
% cost_red_perc_11000NE = 100-(abs(cost_ini_11000NE-cost_fin_11000NE)/abs(cost_ini_11000NE) * 100);
% 
% cost_ini_11BR = cost_lin_ini_BR - cost_rot_ini_BR;
% cost_fin_11BR = cost_lin_fin_11BR - cost_rot_fin_11BR;
% cost_red_perc_11BR = 100-(abs(cost_ini_11BR-cost_fin_11BR)/abs(cost_ini_11BR) * 100);
% 
% cost_fin_noDir = cost_lin_fin_noDir - cost_rot_fin_noDir;
% cost_red_perc_noDir = 100-(abs(cost_lin_ini_NE-cost_fin_noDir)/abs(cost_lin_ini_NE) * 100);
% 
% cost_ini_Straight = cost_lin_ini_Straight - cost_rot_ini_Straight;
% cost_fin_11Straight = cost_lin_fin_11Straight - cost_rot_fin_11Straight;
% cost_red_perc_11Straight = 100-(abs(cost_ini_Straight-cost_fin_11Straight)/abs(cost_ini_Straight) * 100);
% cost_red_perc = [cost_red_perc_11NE cost_red_perc_11000NE cost_red_perc_11BR cost_red_perc_11Straight cost_red_perc_noDir];
% 
% Color = [0,0,1; 0,0,0; 1,0,0; 0.1,0.8,1; 0.7,0.1,0.7; 1,0.5,0.5; 0,1,0; 1.,.5,0.];
% 
% figure(1); clf;
% subplot(1,2,1)
% set(gca, 'FontSize', 15)
% h(1) = plot(cost_lin_tot_NE,'color',Color(1,:),'LineWidth',2);
% hold on
% h(2) = plot(cost_lin_tot_BR,'color',Color(2,:),'LineWidth',2);
% hold on
% h(3) = plot(cost_lin_tot_Straight,'color',Color(3,:),'LineWidth',2);
% hold on
% h(4) = plot(cost_lin_tot_11NE,'color',Color(4,:),'LineWidth',2);
% hold on
% h(5) = plot(cost_lin_tot_11000NE,'color',Color(5,:),'LineWidth',2);
% hold on
% h(6) = plot(cost_lin_tot_11BR,'color',Color(6,:),'LineWidth',2);
% hold on
% h(7) = plot(cost_lin_tot_11Straight,'color',Color(7,:),'LineWidth',2);
% hold on
% h(8) = plot(cost_lin_tot_noDir,'color',Color(8,:),'LineWidth',2);
% hold on
% 
% legend(h,'Neoprene - Initial Shape','Butyl Rubber - Initial Shape','Straight - Initial Shape','Case 1','Case 2','Case 3','Case 4','Case 5','Location','northwest');
% 
% xlabel('Time Step','fontsize',15);
% ylabel('\sigma_{max} (K_{c}^{tr}) [N/m]','fontsize',15);
% ite = size(cost_lin_tot_11000NE,2);
% ite_lin = round(ite/10);
% x = [ite_lin;ite_lin];
% y = [0;18*10^5];
% plot(x,y,'r:','LineWidth',2)
% axis([0 size(Kcart_NE_vec,1) 0 20*10^5])
% 
% figure(1);
% subplot(1,2,2)
% set(gca, 'FontSize', 15)
% h(1) = plot(cost_rot_tot_NE,'color',Color(1,:),'LineWidth',2);
% hold on
% h(2) = plot(cost_rot_tot_BR,'color',Color(2,:),'LineWidth',2);
% hold on
% h(3) = plot(cost_rot_tot_Straight,'color',Color(3,:),'LineWidth',2);
% hold on
% h(4) = plot(cost_rot_tot_11NE,'color',Color(4,:),'LineWidth',2);
% hold on
% h(5) = plot(cost_rot_tot_11000NE,'color',Color(5,:),'LineWidth',2);
% hold on
% h(6) = plot(cost_rot_tot_11BR,'color',Color(6,:),'LineWidth',2);
% hold on
% h(7) = plot(cost_rot_tot_11Straight,'color',Color(7,:),'LineWidth',2);
% hold on
% h(8) = plot(cost_rot_tot_noDir,'color',Color(8,:),'LineWidth',2);
% hold on
% 
% legend(h,'Neoprene - Initial Shape','Butyl Rubber - Initial Shape','Straight - Initial Shape','Case 1','Case 2','Case 3','Case 4','Case 5','Location','northwest');
% xlabel('Time Step','fontsize',15);
% ylabel('\sigma_{min} (K_{c}^{rot}) [N/rad]','fontsize',15);
% axis([0 size(Kcart_NE_vec,1) 0 3000])
% 
% %%% Plot the reduction of the cost function
% figure(2); clf
% 
% set(gca, 'FontSize', 25)
% % Create a vertical bar chart using the bar function
% cost_ini = [100 100 100 100 100];
% h(1) = bar(1:5, cost_ini', 'grouped', 'yellow');
% hold on
% h(2) = bar(1:5, cost_red_perc', 'grouped', 'blue');
% legend(h,'Initial cost function','% of final cost function compared to the initial cost function','Location','north');
% % Set the axis limits
% axis([0 6 0 150])
% 
% % Add title and axis labels
% % title('Reduction of the cost function compared to the initial shape','fontsize',15)
% xlabel('Case','fontsize',25)
% ylabel('%','fontsize',25)
% 
% % Add a legend
% % legend('Measles')



format long
clc
clear all
close all

% addpath ./results/CMAME/CartesianStiffness/ButylRubber
% addpath ./results/CMAME/CartesianStiffness/Neoprene
% addpath ./results/CMAME/CartesianStiffness/NeopreneStraight
% addpath ./results/CMAME/CartesianStiffness/1-1ButylRubber
% addpath ./results/CMAME/CartesianStiffness/1-1Neoprene
% addpath ./results/CMAME/CartesianStiffness/1-1NeoprenenoDir
% addpath ./results/CMAME/CartesianStiffness/1-1NeopreneStraight
% addpath ./results/CMAME/CartesianStiffness/1-1000Neoprene

a = load('K_cart_opt.mat');
K_cart_opt_vec = a.K_cart_vec;
cost_rot_tot_opt = zeros(1,size(K_cart_opt_vec,1));
cost_lin_tot_opt = zeros(1,size(K_cart_opt_vec,1));
a = load('K_cart_ini.mat');
K_cart_ini_vec = a.K_cart_vec;
cost_rot_tot_ini = zeros(1,size(K_cart_ini_vec,1));
cost_lin_tot_ini = zeros(1,size(K_cart_ini_vec,1));
for i=1:size(K_cart_opt_vec,1)
    K_cart_opt = K_cart_opt_vec{i};
    cost_rot_tot_opt(i) = sqrt(min(eig(K_cart_opt(4:6,4:6)'*K_cart_opt(4:6,4:6))));
    cost_lin_tot_opt(i) = sqrt(max(eig(K_cart_opt(1:3,1:3)'*K_cart_opt(1:3,1:3))));
    K_cart_ini = K_cart_ini_vec{i};
    cost_rot_tot_ini(i) = sqrt(min(eig(K_cart_ini(4:6,4:6)'*K_cart_ini(4:6,4:6))));
    cost_lin_tot_ini(i) = sqrt(max(eig(K_cart_ini(1:3,1:3)'*K_cart_ini(1:3,1:3))));     
end
cost_lin_opt = 0;
cost_lin_ini = 0;
for i=1:(size(K_cart_opt_vec,1)/10)
    cost_lin_opt = cost_lin_opt + cost_lin_tot_opt(i);
    cost_lin_ini = cost_lin_ini + cost_lin_tot_ini(i);
end
cost_rot_opt = sum(cost_rot_tot_opt);
cost_rot_ini = sum(cost_rot_tot_ini);
w = 1;
cost_ini = cost_lin_ini - w*cost_rot_ini;
cost_opt = cost_lin_opt - w*cost_rot_opt;


cost_red_perc = 100-(abs(cost_ini-cost_opt)/abs(cost_ini) * 100);

Color = [0,0,1; 0,0,0; 1,0,0; 0.1,0.8,1; 0.7,0.1,0.7; 1,0.5,0.5; 0,1,0; 1.,.5,0.];

figure(1); clf;
subplot(1,2,1)
set(gca, 'FontSize', 15)
h(1) = plot(cost_lin_tot_ini,'color',Color(1,:),'LineWidth',2);
hold on
h(2) = plot(cost_lin_tot_opt,'color',Color(2,:),'LineWidth',2);
hold on

legend(h,'Neoprene - Initial Shape','Optimized');

xlabel('Time Step','fontsize',15);
ylabel('\sigma_{max} (K_{c}^{tr}) [N/m]','fontsize',15);
ite = size(cost_lin_tot_ini,2);
ite_lin = round(ite/10);
x = [ite_lin;ite_lin];
y = [0;18*10^5];
plot(x,y,'r:','LineWidth',2)
axis([0 size(K_cart_ini_vec,1) 0 20*10^5])

figure(1);
subplot(1,2,2)
set(gca, 'FontSize', 15)
h(1) = plot(cost_rot_tot_ini,'color',Color(1,:),'LineWidth',2);
hold on
h(2) = plot(cost_rot_tot_opt,'color',Color(2,:),'LineWidth',2);
hold on

legend(h,'Neoprene - Initial Shape','Optimized');

xlabel('Time Step','fontsize',15);
ylabel('\sigma_{min} (K_{c}^{rot}) [N m/rad]','fontsize',15);
axis([0 size(K_cart_ini_vec,1) 0 3000])

%%% Plot the reduction of the cost function
figure(2); clf

set(gca, 'FontSize', 25)
% Create a vertical bar chart using the bar function
cost_ini = [100 100 100 100 100];
h(1) = bar(1:5, cost_ini', 'grouped', 'yellow');
hold on
h(2) = bar(1:5, cost_red_perc', 'grouped', 'blue');
legend(h,'Initial cost function','% of final cost function compared to the initial cost function','Location','north');
% Set the axis limits
axis([0 6 0 150])

% Add title and axis labels
% title('Reduction of the cost function compared to the initial shape','fontsize',15)
xlabel('Case','fontsize',25)
ylabel('%','fontsize',25)

% Add a legend
% legend('Measles')

