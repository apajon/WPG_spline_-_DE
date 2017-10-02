function p = reshape_p_vec(p_vec,spl,row_spli,col_spli)

assert(row_spli==col_spli,'Not implemented for row_spli1~=col_spli')

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
end