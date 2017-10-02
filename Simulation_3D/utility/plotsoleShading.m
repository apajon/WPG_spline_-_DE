function fig = plotsoleShading(number,connecext,coor,stressVM,az,el)
fig = figure(number);
%set(fig,'PlotBoxAspectRatioMode','Manual')
% pbaspect([0.5,0.5,0.5])
clf();
a = patch('faces',connecext,'vertices',coor,'FaceVertexCData',stressVM(:,1),'FaceColor','yellow','EdgeColor','none');
% set(a,'DataAspectRatio',[1 1 1],...
%         'PlotBoxAspectRatio',[1 1 1],'ZLim',[-0.6 0.6])
%camlight left;
%shading interp
% zoom off

camva('manual')
camva(10)
% set(fig,'PlotBoxAspectRatioMode','Manual')
% pbaspect([1.0,1.0,1.0])
% lighting phong
% lighting flat
lighting gouraud
h1 = light('Position',[-0.1 0 0.1],'Style','infinite');
% h1.color = [1 1 0];
h2 = light('Position',[0.3 0 -0.1],'Style','infinite');
% h2.color = [1 1 0];
% camlight(-60,30)
% camlight(-120,30)
% shading interp
% h = camlight('left');
% for i = 1:20;
%    camorbit(10,0)
%    camlight(h,'left')
%    pause(.1)
% end
% grid on
% axis off
axis([-0.1  0.3 -0.24  0.24  -0.2  0.15])
view([az,el])

xlabel('x[m]') % x-axis label
ylabel('y[m]') % y-axis label
zlabel('z[m]') % z-axis label

hold on;
% pointA = [-10,10,0];
% pointB = [-10,-10,0];
% pointC = [10,-10,0];pointD = [10,10,0];
% points = [pointA' pointB' pointC' pointD']; % using the data given in the question
% fill3(points(1,:),points(2,:),points(3,:),'r')
% alpha(1)
% hold on;
%view(90,0)

% [X,Y] = meshgrid(-0.3:.05:0.3);
% R = sqrt(X.^2 + Y.^2) + eps;
% Z = sin(0)./R;
% % patch(surf2patch(X,Y,Z,Z),'FaceColor','g','EdgeColor','flat',...
% %       'Marker','o','MarkerFaceColor','flat');
% colorPlane = zeros(size(X,1)*size(X,1),3);
% colorPlane(:,1) = 0.9;
% colorPlane(:,2) = 0.9;
% colorPlane(:,3) = 0.9;
% patch(surf2patch(X,Y,Z,Z),'FaceVertexCData',colorPlane);
% shading faceted;
% hold on;
%view(3)
%In case of contact
axis image
% daspect manual
end
