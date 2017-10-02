function sole = fitting_parallelepiped(sole,l,L,e,Point)
    %%% fit inverse data to the original 3d parallelepiped model
    nodesSurf=sole.coor(sole.nodesSurf,:);
    z_max = max(nodesSurf(:,3));
    z_min = min(nodesSurf(:,3));
    sole.coor(sole.nodesDirichlet,3)=z_max;  

    %B = setdiff(sole.nodesFreeSurf,sole.nodesSurf(I_surf_bottom));
    B = sole.nodesFreeSurf;
    figure(4); clf; plot3(sole.coor(B,1),sole.coor(B,2),sole.coor(B,3),'r+')
    plane_yz_back = createPlane(Point(1,:),Point(2,:),Point(3,:));
    n_yz_back = planeNormal(plane_yz_back);
    plane_yz_front = createPlane(Point(5,:),Point(6,:),Point(7,:));
    n_yz_front = planeNormal(plane_yz_front);    
    plane_xz_right = createPlane(Point(2,:),Point(3,:),Point(6,:));
    n_xz_right = planeNormal(plane_xz_right);
    plane_xz_left = createPlane(Point(1,:),Point(4,:),Point(5,:));
    n_xz_left = planeNormal(plane_xz_left);    
    plane_xy_bottom = createPlane(Point(1,:),Point(2,:),Point(6,:));
    n_xy_bottom = planeNormal(plane_xy_bottom);
    
%     plotsole(6,sole.elements_surf,sole.coor,stressVM0,[],[],[],[],-37.5,30);  
%     hold on
%     figure(6); clf; drawPlane3d(plane_xz_right)
    for i=1:size(B,1)
        %%% yz back
        [coor_yz_back,check_yz_back]=plane_line_intersect(n_yz_back,Point(1,:),[L/2,0,e],sole.coor(B(i),:));
        distancePoint_yz_back = sqrt((sole.coor(B(i),1)-coor_yz_back(1))^2+(sole.coor(B(i),2)-coor_yz_back(2))^2+(sole.coor(B(i),3)-coor_yz_back(3))^2);
        
        %%% yz front
        [coor_yz_front,check_yz_front]=plane_line_intersect(n_yz_front,Point(5,:),[L/2,0,e],sole.coor(B(i),:));
        distancePoint_yz_front = sqrt((sole.coor(B(i),1)-coor_yz_front(1))^2+(sole.coor(B(i),2)-coor_yz_front(2))^2+(sole.coor(B(i),3)-coor_yz_front(3))^2);
        
         %%% xz right
        [coor_xz_right,check_xz_right]=plane_line_intersect(n_xz_right,Point(6,:),[L/2,0,e],sole.coor(B(i),:));
        distancePoint_xz_right = sqrt((sole.coor(B(i),1)-coor_xz_right(1))^2+(sole.coor(B(i),2)-coor_xz_right(2))^2+(sole.coor(B(i),3)-coor_xz_right(3))^2);    
        
        %%% xz left
        [coor_xz_left,check_xz_left]=plane_line_intersect(n_xz_left,Point(5,:),[L/2,0,e],sole.coor(B(i),:));
        distancePoint_xz_left = sqrt((sole.coor(B(i),1)-coor_xz_left(1))^2+(sole.coor(B(i),2)-coor_xz_left(2))^2+(sole.coor(B(i),3)-coor_xz_left(3))^2);     
        
        %%% xy bottom
        [coor_xy_bottom,check_xy_bottom]=plane_line_intersect(n_xy_bottom,Point(6,:),[L/2,0,e],sole.coor(B(i),:));
        distancePoint_xy_bottom = sqrt((sole.coor(B(i),1)-coor_xy_bottom(1))^2+(sole.coor(B(i),2)-coor_xy_bottom(2))^2+(sole.coor(B(i),3)-coor_xy_bottom(3))^2);
      
        
        coor_up = [coor_yz_back;coor_yz_front;coor_xz_right;coor_xz_left;coor_xy_bottom];
        distance = [distancePoint_yz_back,distancePoint_yz_front,distancePoint_xz_right,distancePoint_xz_left,distancePoint_xy_bottom];
        
        %%% find the closest plane and project the coordinate in this plane
        [min_dist,I_dist_min]=min(distance);
        %if (distancePoint_xz_left < 5e-2)
        sole.coor(B(i),:) = coor_up(I_dist_min,:);
        %end
%         hold on
%         figure(6); plot3(sole.coor(B(i),1),sole.coor(B(i),2),sole.coor(B(i),3),'g+')
%         hold on         
        
    end
%     figure(6); plot3(sole.coor(B,1),sole.coor(B,2),sole.coor(B,3),'r+')

end