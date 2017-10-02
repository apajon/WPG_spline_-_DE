function convert_to_gmesh_file(sole,p_ini_v,spl,param_sopt,fname)
    sole.coor = deformation_moveDiri(sole,p_ini_v,spl,param_sopt.move_dirichlet); % Moving Dirichlet

    coorNodesnew(:,1) = 1:sole.nTot;
    coorNodesnew(:,2) = sole.coor(:,1);
    coorNodesnew(:,3) = sole.coor(:,2);
    coorNodesnew(:,4) = sole.coor(:,3);
    coorNodesnew = coorNodesnew';
    %%% Define surfaace elements
    eleSurf(:,1) = 1:size(sole.elements_surf,1);
    eleSurf(:,2) = sole.elements_surf(:,1);
    eleSurf(:,3) = sole.elements_surf(:,2);
    eleSurf(:,4) = sole.elements_surf(:,3);
    eleSurf = eleSurf';   
    nodeTriangleType = zeros(size(sole.elements_surf,1),1)+2; % 2 = 3-node triangle type
    physicalSurf = zeros(size(sole.elements_surf,1),1)+10;
    elementarySurf = zeros(size(sole.elements_surf,1),1)+2;
    %%% Define Dirichlet surface for gmesh
    elementarySurf1 = zeros(size(sole.elements_surf,1),1);
    LIA = ismember(sole.elements_surf,sole.nodesDirichlet);
    for i=1:size(sole.elements_surf,1)
        elementarySurf1(i) = 6;
        if LIA(i,1)==1 && LIA(i,2)==1 && LIA(i,3)==1
            elementarySurf1(i) = 1;
        end
    end
    %%% Define volume elements
    formatSurface = [eleSurf(1,:);nodeTriangleType';elementarySurf';physicalSurf';elementarySurf1';eleSurf(2:4,:)];
    eleVol(:,1) = (1:size(sole.elements_vol,1))+size(sole.elements_surf,1);
    eleVol(:,2) = sole.elements_vol(:,1);
    eleVol(:,3) = sole.elements_vol(:,2);
    eleVol(:,4) = sole.elements_vol(:,3);
    eleVol(:,5) = sole.elements_vol(:,4);
    eleVol = eleVol';  
    nodeTetrahedronType = zeros(size(sole.elements_vol,1),1)+4; % 4 = 4-node tetrahedron.
    elementaryVol = zeros(size(sole.elements_vol,1),1)+2;
    physicalVol = zeros(size(sole.elements_vol,1),1)+100;
    formatVolume = [eleVol(1,:);nodeTetrahedronType';elementaryVol';physicalVol';elementaryVol';eleVol(2:5,:)];
    
    fid = fopen(fname,'w');
    fprintf(fid,'$MeshFormat \n2.2 0 8\n$EndMeshFormat\n$Nodes\n%d\n',sole.nTot);
    fprintf(fid,'%d %f %f %f\n',coorNodesnew);
    fprintf(fid,'$EndNodes\n');
    fprintf(fid,'$Elements\n');
    fprintf(fid,'%d\n', sole.nEle);
    fprintf(fid,'%d %d %d %d %d %d %d %d\n',formatSurface);
    fprintf(fid,'%d %d %d %d %d %d %d %d %d\n',formatVolume);
    fprintf(fid,'$EndElements\n');
    fclose(fid);
end