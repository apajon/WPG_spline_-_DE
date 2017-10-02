function fig = plotsoleShading(number,connecext,coor,stressVM,az,el)
fig = figure(number);
clf();
a = patch('faces',connecext,'vertices',coor,'FaceVertexCData',stressVM(:,1),'FaceColor','yellow','EdgeColor','none');

camva('manual')
camva(10)
lighting gouraud
h1 = light('Position',[-0.1 0 0.1],'Style','infinite');
h2 = light('Position',[0.3 0 -0.1],'Style','infinite');

axis([-0.1  0.3 -0.24  0.24  -0.2  0.15])
view([az,el])

set(gca, 'FontSize', 30)
xlabel('x[m]','fontsize',30);
ylabel('y[m]','fontsize',30);
zlabel('z[m]','fontsize',30);

hold on;
axis image

end
