function coornew = deformation_moveDiri_inverse(sole,p_vec,spl,move_dirichlet)
s_coor = size(sole.coor);
contSurf = 1;
if move_dirichlet == 1
    p_mid = reshape(p_vec(1:(spl.l_spli*spl.l_spli)),spl.l_spli,spl.l_spli);
    iniv1 = (spl.l_spli*spl.l_spli)+1;
    ptmp = p_vec(iniv1:(iniv1+spl.l_spli-2-1))';

    pfirstpart = p_vec(end-1);
    plastpart = p_vec(end);
    p1 = [pfirstpart pfirstpart ptmp plastpart plastpart];
    iniv2 = iniv1+spl.l_spli-2;
    ptmp2 = p_vec(iniv2:(iniv2+spl.l_spli-2-1))';
    p2 = [pfirstpart pfirstpart ptmp2 plastpart plastpart];
    p3 = zeros(size(p_mid,1),1)+pfirstpart; 
    p4 = zeros(size(p_mid,1),1)+plastpart;
    p = [p1; p3 p_mid p4; p2];
    coornew = zeros(s_coor);
    coor_tmp = zeros(1,3);
    for i=1:length(spl.r)
        s1=sum(sum(p(spl.idx_u{i}, spl.idx_v{i}).*spl.basisValues{i}));
        [coor_tmp(1,1),coor_tmp(1,2),coor_tmp(1,3)] = sph2cart(spl.u{i},spl.v{i},(1/s1)*spl.r{i});
        coornew(i,:) = [coor_tmp(1,3),coor_tmp(1,2),-coor_tmp(1,1)];
    end
else
    p_mid = reshape(p_vec(1:(spl.l_spli*spl.l_spli)),spl.l_spli,spl.l_spli);
    p = [ones(1,size(p_mid,1)+2); ones(size(p_mid,1),1) p_mid ones(size(p_mid,1),1); ones(1,size(p_mid,1)+2)];
    coornew = zeros(s_coor);
    coor_tmp = zeros(1,3);
    for i=1:length(spl.r)
        s1=sum(sum(p(spl.idx_u{i}, spl.idx_v{i}).*spl.basisValues{i}));
        [coor_tmp(1,1),coor_tmp(1,2),coor_tmp(1,3)] = sph2cart(spl.u{i},spl.v{i},(1/s1)*spl.r{i});
        coornew(i,:) = [coor_tmp(1,3),coor_tmp(1,2),-coor_tmp(1,1)];    
    end
end
%%% translate because spl is in coornew-sole.trasl
coornew(:,1)= coornew(:,1) + sole.trasl(1);
coornew(:,2) = coornew(:,2) + sole.trasl(2);
coornew(:,3) = coornew(:,3) + sole.trasl(3);
end