function sole = remesh_isomesh_soleini(sole,k,reduce,remesh,l,L,e)
contSurf=1;
for i=1:sole.nTot
    if contSurf<=length(sole.nodesSurf)
        if i==sole.nodesSurf(contSurf)
            I1 = find(sole.elements_surf(:,1)==i);
            sole.elements_surf(I1,1)=contSurf;
            I2 = find(sole.elements_surf(:,2)==i);
            sole.elements_surf(I2,2)=contSurf;
            I3 = find(sole.elements_surf(:,3)==i);
            sole.elements_surf(I3,3)=contSurf;
            contSurf = contSurf + 1;
        end
    end
end
a = sole.coor(sole.nodesSurf,:);
b = sole.elements_surf;
if reduce == 1
   % reduce the surface node numbers to (p*100)%
   p = 0.5;
   [a,b] = meshresample(a,b,p);  % down sample the sphere mesh
end
% figure; plotmesh(a,b);
if remesh == 1
    [coor1,eleVol1,face1] = surf2mesh(a,b,[0 0 0],[0. 0. 0.],k,1,[0 0 0],[],1);
else
    coor1 = sole.coor;
    eleVol1 = sole.elements_vol;
    face1 = sole.elements_surf;
end

sole.coor = [];
sole.coor = coor1;
sole.elements_vol = [];
sole.elements_vol = [eleVol1(:,1) eleVol1(:,2) eleVol1(:,3) eleVol1(:,4)];
sole.elements_surf = [];
sole.elements_surf = [face1(:,1) face1(:,2) face1(:,3)];

nodes = (1:1:size(sole.coor,1))';
sole.lx_foot = L;
sole.ly_foot = l;
sole.lz_foot = e;
sole.nTot = length(nodes);
sole.trasl = [sole.lx_foot/2,0,sole.lz_foot];
% Dirichlet nodes
sole.nodesDirichlet = [];
sole.nodesDirichlet = find(sole.coor(:,3)==sole.lz_foot); 
sole.nDirichlet = length(sole.nodesDirichlet);
sole.nodesSurf = [];
sole.nodesSurf = unique(sole.elements_surf);
% Free surface nodes
sole.nodesFreeSurf = [];
sole.nodesFreeSurf = setdiff(sole.nodesSurf,sole.nodesDirichlet);
sole.nFreeSurf = length(sole.nodesFreeSurf);
% Nodes of volume - Internal nodes
sole.nodesInt = [];
sole.nodesInt = nodes;
sole.nodesInt([sole.nodesFreeSurf; sole.nodesDirichlet])=[];
sole.nInt = length(sole.nodesInt);
% Nodes free - Free surface nodes + Internal nodes
sole.nodesFree = [];
sole.nodesFree = nodes;
sole.nodesFree(sole.nodesDirichlet)=[];
sole.nFree = length(sole.nodesFree);
% Compute dof matrix
cont = 0;
sole.dof = zeros(size(sole.coor));
contDir = 1;
for i=1:sole.nTot
    if contDir <= sole.nDirichlet
        if i~=sole.nodesDirichlet(contDir)
            for j=1:3
                cont = cont + 1;
                sole.dof(i,j) = cont;
            end
        else
            contDir = contDir + 1;
        end
    else
         for j=1:3
            cont = cont + 1;
            sole.dof(i,j) = cont;
        end                   
    end
end
% Nodes in x;y;z;
sole.nodesFree3 = [];
sole.nodesFree3(1:3:3*sole.nFree-2) = 3*sole.nodesFree-2;
sole.nodesFree3(2:3:3*sole.nFree-1) = 3*sole.nodesFree-1;
sole.nodesFree3(3:3:3*sole.nFree-0) = 3*sole.nodesFree-0; 
sole.nodesFreeSurf3 = [];
sole.nodesFreeSurf3(1:3:3*sole.nFreeSurf-2) = 3*sole.nodesFreeSurf-2;
sole.nodesFreeSurf3(2:3:3*sole.nFreeSurf-1) = 3*sole.nodesFreeSurf-1;
sole.nodesFreeSurf3(3:3:3*sole.nFreeSurf-0) = 3*sole.nodesFreeSurf-0;
sole.nodesInt3 = [];
sole.nodesInt3(1:3:3*sole.nInt-2) = 3*sole.nodesInt-2;
sole.nodesInt3(2:3:3*sole.nInt-1) = 3*sole.nodesInt-1;
sole.nodesInt3(3:3:3*sole.nInt-0) = 3*sole.nodesInt-0;
sole.nodesDirichlet3 = [];
sole.nodesDirichlet3(1:3:3*sole.nDirichlet-2) = 3*sole.nodesDirichlet-2;
sole.nodesDirichlet3(2:3:3*sole.nDirichlet-1) = 3*sole.nodesDirichlet-1;
sole.nodesDirichlet3(3:3:3*sole.nDirichlet-0) = 3*sole.nodesDirichlet-0;
