function [fig1,subplot1,subplot2] = drawSimulation(sole,Pg,stressVM,Pc,Fc,Z,Ftot,Zdes,Fdes,i,w,move_dirichlet,az,el)
% figure(2)
% subplot(2,2,4)
% hold on
% plot(Zdes(1,i),Zdes(2,i),'o','LineWidth',25)
% hold off
% 
% figure(2)
% subplot(2,2,3)
% hold on
% plot(i,Fdes(1,i),'ko','LineWidth',25)
% plot(i,Fdes(2,i),'ro','LineWidth',25)
% plot(i,Fdes(3,i),'bo','LineWidth',25)
% 
% hold off 
% fig1 = figure(2);
% % clf();
% subplot1 = subplot(2,2,1); 
% 
% patch('faces',sole.elements_surf,'vertices',Pg,'FaceVertexCData',stressVM,'FaceColor','interp','CDataMapping','scaled');
% % grid on
% axis([-0.1  0.3 -0.24  0.24  -0.2  0.15])
% %axis([-10  30 -24  24  -20  15])
% view([-37.5,30])
% 
% set(gca, 'FontSize', 25)
% xlabel('x[m]','fontsize',25);
% ylabel('y[m]','fontsize',25);
% zlabel('z[m]','fontsize',25);
% 
% hold on;
% %In case of contact
% if ~isempty(Fc)%----------------------------------------------------
%     alpha(1)
%     % Plot the contact forces
%     quiver3(Pc(:,1),Pc(:,2),Pc(:,3),-Fc(:,1)*0.005,-Fc(:,2)*0.005,-Fc(:,3)*0.005,0,'b')
% 
%     %-----------------------calcul ZMP position--------------------------------
%     %Fsol = sum(Fc,1);
%     Mo = sum(cross(Pc',Fc'),2);
%     Xzmp = Z(1);
%     Yzmp = Z(2);
%     Zzmp = 0;
% 
%     % Track the Zmp (green line)
%     %quiver3(Xzmp,Yzmp,Zzmp,-Ftot(1)*0.0003,-Ftot(2)*0.0003,-Ftot(3)*0.0003,0)
%     
%     drawnow();
% end
% 
% subplot2 = subplot(2,2,2);
% patch('faces',sole.elements_surf,'vertices',Pg,'FaceVertexCData',stressVM,'FaceColor','interp','CDataMapping','scaled');
% % grid on
% axis([-0.1  0.3 -0.24  0.24  -0.2  0.15])
% %axis([-10  30 -24  24  -20  15])
% view([-127.5,30])
% 
% set(gca, 'FontSize', 25)
% xlabel('x[m]','fontsize',25);
% ylabel('y[m]','fontsize',25);
% zlabel('z[m]','fontsize',25);
% 
% hold on;
% %In case of contact
% if ~isempty(Fc)%----------------------------------------------------
%     alpha(1)
%     % Plot the contact forces
%     quiver3(Pc(:,1),Pc(:,2),Pc(:,3),-Fc(:,1)*0.005,-Fc(:,2)*0.005,-Fc(:,3)*0.005,0,'b')
% 
%     %-----------------------calcul ZMP position--------------------------------
%     %Fsol = sum(Fc,1);
%     Mo = sum(cross(Pc',Fc'),2);
%     Xzmp = Z(1);
%     Yzmp = Z(2);
%     Zzmp = 0;
% 
%     % Track the Zmp (green line)
%     %quiver3(Xzmp,Yzmp,Zzmp,-Ftot(1)*0.0003,-Ftot(2)*0.0003,-Ftot(3)*0.0003,0)
%     
%     drawnow();
% end
% 
% set(gcf,'NextPlot','add');
% axes;
% if w==0
%     t1 = '\fontsize{55}Simulation - Initial sole shape:        Material=';
% else
%     t1 = '\fontsize{55}Simulation - Optimized sole shape:        Material=';
% end
% if sole.E==1000000
%     t2 = 'Neoprene';
% else
%     t2 = 'Butyl Rubber';
% end
% if w==0
%     t3='';
%     t4='';
% else    
%     if w == 1
%         t3 = ', Cost: w=1, ';
%     else
%         t3 = ', Cost: w=1000, ';
%     end
%     if move_dirichlet==1
%         t4 = 'Dirichlet=Yes';
%     else
%         t4 = 'Dirichlet=No';
%     end
% end
% t = strcat(t1,t2,t3,t4); 
% h = title(t);
% set(gca,'Visible','off');
% set(h,'Visible','on');
subplot1 = [];
subplot2 = [];
figure(2)
% % subplot(2,2,4)
% hold on
% plot(Zdes(1,i),Zdes(2,i),'o','LineWidth',14)
% hold off

% figure(2)
% subplot(2,2,3)
% hold on
% plot(i,Fdes(1,i),'ko','LineWidth',25)
% plot(i,Fdes(2,i),'ro','LineWidth',25)
% plot(i,Fdes(3,i),'bo','LineWidth',25)
% 
% hold off 
fig1 = figure(2);
clf();
% subplot1 = subplot(2,2,1); 
% 
patch('faces',sole.elements_surf,'vertices',Pg,'FaceVertexCData',stressVM,'FaceColor','interp','CDataMapping','scaled');
% grid on
axis([-0.1  0.3 -0.24  0.24  -0.2  0.15])
%axis([-10  30 -24  24  -20  15])
view([az,el])

set(gca, 'FontSize', 16)
xlabel('x[m]','fontsize',16);
ylabel('y[m]','fontsize',16);
zlabel('z[m]','fontsize',16);

hold on;
%In case of contact
if ~isempty(Fc)%----------------------------------------------------
    alpha(1)
    % Plot the contact forces
    quiver3(Pc(:,1),Pc(:,2),Pc(:,3),-Fc(:,1)*0.005,-Fc(:,2)*0.005,-Fc(:,3)*0.005,0,'b')

    %-----------------------calcul ZMP position--------------------------------
    %Fsol = sum(Fc,1);
    Mo = sum(cross(Pc',Fc'),2);
    Xzmp = Z(1);
    Yzmp = Z(2);
    Zzmp = 0;

    % Track the Zmp (green line)
    %quiver3(Xzmp,Yzmp,Zzmp,-Ftot(1)*0.0003,-Ftot(2)*0.0003,-Ftot(3)*0.0003,0)
    
    drawnow();
end
title('\fontsize{16} Neoprene, Energy-minimizing trajectory, w=1, Fixed Dirichlet nodes')
% subplot2 = subplot(2,2,2);
% patch('faces',sole.elements_surf,'vertices',Pg,'FaceVertexCData',stressVM,'FaceColor','interp','CDataMapping','scaled');
% % grid on
% axis([-0.1  0.3 -0.24  0.24  -0.2  0.15])
% %axis([-10  30 -24  24  -20  15])
% view([-127.5,30])
% 
% set(gca, 'FontSize', 25)
% xlabel('x[m]','fontsize',25);
% ylabel('y[m]','fontsize',25);
% zlabel('z[m]','fontsize',25);
% 
% hold on;
% %In case of contact
% if ~isempty(Fc)%----------------------------------------------------
%     alpha(1)
%     % Plot the contact forces
%     quiver3(Pc(:,1),Pc(:,2),Pc(:,3),-Fc(:,1)*0.005,-Fc(:,2)*0.005,-Fc(:,3)*0.005,0,'b')
% 
%     %-----------------------calcul ZMP position--------------------------------
%     %Fsol = sum(Fc,1);
%     Mo = sum(cross(Pc',Fc'),2);
%     Xzmp = Z(1);
%     Yzmp = Z(2);
%     Zzmp = 0;
% 
%     % Track the Zmp (green line)
%     %quiver3(Xzmp,Yzmp,Zzmp,-Ftot(1)*0.0003,-Ftot(2)*0.0003,-Ftot(3)*0.0003,0)
%     
%     drawnow();
% end
% 
% set(gcf,'NextPlot','add');
% axes;
% if w==0
%     t1 = '\fontsize{55}Simulation - Initial sole shape:        Material=';
% else
%     t1 = '\fontsize{55}Simulation - Optimized sole shape:        Material=';
% end
% if sole.E==1000000
%     t2 = 'Neoprene';
% else
%     t2 = 'Butyl Rubber';
% end
% if w==0
%     t3='';
%     t4='';
% else    
%     if w == 1
%         t3 = ', Cost: w=1, ';
%     else
%         t3 = ', Cost: w=1000, ';
%     end
%     if move_dirichlet==1
%         t4 = 'Dirichlet=Yes';
%     else
%         t4 = 'Dirichlet=No';
%     end
% end
% t = strcat(t1,t2,t3,t4); 
% h = title(t);
% set(gca,'Visible','off');
% set(h,'Visible','on');



end