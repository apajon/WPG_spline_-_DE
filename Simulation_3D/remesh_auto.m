function [sole,spl] = remesh_auto(sole,p_ini_v,spl,move_dirichlet,l,L,e)
    sole.coor = deformation_moveDiri(sole,p_ini_v,spl,move_dirichlet);    
    stressVM0 = zeros(sole.nTot,1);
    plotsole(1,sole.elements_surf,sole.coor,stressVM0,[],[],[],[],-37.5,30);     
    remesh = 1;
    reduce = 0;
    k = 0.9;
    sole = remesh_isomesh_soleini(sole,k,reduce,remesh,l,L,e);    
    stressVM0 = zeros(sole.nTot,1);
    plotsole(5,sole.elements_surf,sole.coor,stressVM0,[],[],[],[],-37.5,30);      
 
    %%% B-spline remesh coordinates
    spl.update(sole);
    sole.coor = deformation_moveDiri_inverse(sole,p_ini_v,spl,move_dirichlet);
  
    spl.update(sole);
    
    stressVM0 = zeros(sole.nTot,1);
    plotsole(1,sole.elements_surf,sole.coor,stressVM0,[],[],[],[],-37.5,30);
end