function draw3dview(sole,p_ini_v,spl,param_sopt,s1,move_dirichlet,w)
sole.coor = deformation_moveDiri(sole,p_ini_v,spl,param_sopt.move_dirichlet); % Moving Dirichlet

stressVM0 = zeros(sole.nTot,1);
contImg = 0;
t = '\fontsize{16} Neoprene, Energy-minimizing trajectory, w=1000, Free Dirichlet nodes';
%t = '\fontsize{16} Unoptimized sole shape';
for i=-40:1:180
    contImg = contImg + 1;
    fig1 = plotsoleShading(1,sole.elements_surf,sole.coor,stressVM0,i,15);
    set(gcf,'NextPlot','add');
    axes;
    %t1 = '\fontsize{55}3D View - Optimized sole shape:        Material=';
%     t1 = '\fontsize{20}3D View - Unoptimized sole shape:';
%     if sole.E==1000000
%         t2 = 'Neoprene';
%     else
%         t2 = 'Butyl Rubber';
%     end
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
%     t = strcat(t1,t2,t3,t4); 
    h = title(t);
    set(gca,'Visible','off');
    set(h,'Visible','on');    
    s2 = num2str(contImg);
    s = strcat(s1,s2); 
%     set(gcf,'PaperUnits','inches','PaperPosition',[0 0 40 30])
    print(fig1,s,'-dpng')
end
for i=15:-1:-90
    contImg = contImg + 1;
    fig1 = plotsoleShading(1,sole.elements_surf,sole.coor,stressVM0,180,i);
    set(gcf,'NextPlot','add');
    axes;
    %t1 = '\fontsize{55}3D View - Optimized sole shape:        Material=';
    %t1 = '\fontsize{55}3D View - Unoptimized sole shape:';
%     if sole.E==1000000
%         t2 = 'Neoprene';
%     else
%         t2 = 'Butyl Rubber';
%     end
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
%     t = strcat(t1,t2,t3,t4); 
    h = title(t);
    set(gca,'Visible','off');
    set(h,'Visible','on');    
    s2 = num2str(contImg);
    s = strcat(s1,s2); 
%     set(gcf,'PaperUnits','inches','PaperPosition',[0 0 40 30])
    print(fig1,s,'-dpng')
end