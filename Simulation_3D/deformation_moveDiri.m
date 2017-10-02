function [coornew,G_spline_poly] = deformation_moveDiri(sole,p_vec,spl,move_dirichlet)
G_spline_poly = zeros(3*sole.nTot,length(p_vec));
if move_dirichlet == 1
    p = reshape_p_vec(p_vec,spl,spl.l_spli+2,spl.l_spli+2);
    coornew = zeros(length(spl.r),3);
    coor_tmp = zeros(1,3);
    for i=1:length(spl.r)
        s1=sum(sum(p(spl.idx_u{i}, spl.idx_v{i}).*spl.basisValues{i}));
        [coor_tmp(1,1),coor_tmp(1,2),coor_tmp(1,3)] = sph2cart(spl.u{i},spl.v{i},s1*spl.r{i});
        [coor_ini_tmp(1,1),coor_ini_tmp(1,2),coor_ini_tmp(1,3)] = sph2cart(spl.u{i},spl.v{i},spl.r{i});
        coornew(i,:) = [coor_tmp(1,3),coor_tmp(1,2),-coor_tmp(1,1)]; 
        %%% Gradient
        idxu = spl.idx_u{i};
        idxv = spl.idx_v{i};
        BV = spl.basisValues{i};
        for j=1:length(idxu)
            for k=1:length(idxv)
                %ind = find_ind_p(idxu(j),idxv(k),length(p_vec),size(p,1));
                %G_spline_poly((3*i-2):(3*i),ind) = G_spline_poly((3*i-2):(3*i),ind) + BV(j,k);
                p_grad = reshape_p_vec([1:length(p_vec)]',spl,spl.l_spli+2,spl.l_spli+2);
                ind = p_grad(idxu(j),idxv(k));
                G_spline_poly((3*i-2):(3*i),ind) = G_spline_poly((3*i-2):(3*i),ind) + (BV(j,k) * (sole.coor(i,:)-sole.trasl)');
            end
        end
    end
else
    p_mid = reshape(p_vec(1:(spl.l_spli*spl.l_spli)),spl.l_spli,spl.l_spli);
    p_mid(10,1:3)=0.85;
    p = [ones(1,size(p_mid,1)+2); ones(size(p_mid,1),1) p_mid ones(size(p_mid,1),1); ones(1,size(p_mid,1)+2)];
    coornew = zeros(length(spl.r),3);
    coor_tmp = zeros(1,3);
    for i=1:length(spl.r)
        s1=sum(sum(p(spl.idx_u{i}, spl.idx_v{i}).*spl.basisValues{i}));
        [coor_tmp(1,1),coor_tmp(1,2),coor_tmp(1,3)] = sph2cart(spl.u{i},spl.v{i},s1*spl.r{i});
        coornew(i,:) = [coor_tmp(1,3),coor_tmp(1,2),-coor_tmp(1,1)];    
        %%% Gradient
        %%% Not Implemented
%         idxu = spl.idx_u{i};
%         idxv = spl.idx_v{i};
%         BV = spl.basisValues{i};
%         for j=1:length(idxu)
%             for k=1:length(idxv)
%                 %ind = find_ind_p(idxu(j),idxv(k),length(p_vec),size(p,1));
%                 %G_spline_poly((3*i-2):(3*i),ind) = G_spline_poly((3*i-2):(3*i),ind) + BV(j,k);
%                 p_grad = reshape_p_vec([1:length(p_vec)]',spl,spl.l_spli+2,spl.l_spli+2);
%                 ind = p_grad(idxu(j),idxv(k));
%                 G_spline_poly((3*i-2):(3*i),ind) = G_spline_poly((3*i-2):(3*i),ind) + (BV(j,k) * (sole.coor(i,:)-sole.trasl)');
%             end
%         end    
    end   
end
%%% translate because spl is in coornew-sole.trasl
coornew(:,1)= coornew(:,1) + sole.trasl(1);
coornew(:,2) = coornew(:,2) + sole.trasl(2);
coornew(:,3) = coornew(:,3) + sole.trasl(3);
% for i=1:length(spl.r)
%     G_spline_poly((3*i-2):(3*i),ind) = G_spline_poly((3*i-2):(3*i),ind) .* coorini(i,:)';
% end
end