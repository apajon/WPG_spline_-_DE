function dvol_dPfree = dvoldPfree(V,dPFree_dpi_v)
x1=V(1,1); dx1_pi=dPFree_dpi_v(1,1);
y1=V(1,2); dy1_pi=dPFree_dpi_v(1,2);
z1=V(1,3); dz1_pi=dPFree_dpi_v(1,3);
x2=V(2,1); dx2_pi=dPFree_dpi_v(2,1);
y2=V(2,2); dy2_pi=dPFree_dpi_v(2,2);
z2=V(2,3); dz2_pi=dPFree_dpi_v(2,3);
x3=V(3,1); dx3_pi=dPFree_dpi_v(3,1); 
y3=V(3,2); dy3_pi=dPFree_dpi_v(3,2);
z3=V(3,3); dz3_pi=dPFree_dpi_v(3,3);
x4=V(4,1); dx4_pi=dPFree_dpi_v(4,1);
y4=V(4,2); dy4_pi=dPFree_dpi_v(4,2);
z4=V(4,3); dz4_pi=dPFree_dpi_v(4,3);
vol_dx1 = (y3*z2)/6 - (y2*z3)/6 + (y2*z4)/6 - (y4*z2)/6 - (y3*z4)/6 + (y4*z3)/6;
dvolx1_dp_i = vol_dx1 * dx1_pi;
vol_dy1 = (x2*z3)/6 - (x3*z2)/6 - (x2*z4)/6 + (x4*z2)/6 + (x3*z4)/6 - (x4*z3)/6;
dvoly1_dp_i = vol_dy1 * dy1_pi;
vol_dz1 = (x3*y2)/6 - (x2*y3)/6 + (x2*y4)/6 - (x4*y2)/6 - (x3*y4)/6 + (x4*y3)/6;
dvolz1_dp_i = vol_dz1 * dz1_pi;

vol_dx2 = (y1*z3)/6 - (y3*z1)/6 - (y1*z4)/6 + (y4*z1)/6 + (y3*z4)/6 - (y4*z3)/6;
dvolx2_dp_i = vol_dx2 * dx2_pi;
vol_dy2 = (x3*z1)/6 - (x1*z3)/6 + (x1*z4)/6 - (x4*z1)/6 - (x3*z4)/6 + (x4*z3)/6;
dvoly2_dp_i = vol_dy2 * dy2_pi;
vol_dz2 = (x1*y3)/6 - (x3*y1)/6 - (x1*y4)/6 + (x4*y1)/6 + (x3*y4)/6 - (x4*y3)/6;
dvolz2_dp_i = vol_dz2 * dz2_pi;

vol_dx3 = (y2*z1)/6 - (y1*z2)/6 + (y1*z4)/6 - (y4*z1)/6 - (y2*z4)/6 + (y4*z2)/6;
dvolx3_dp_i = vol_dx3 * dx3_pi;
vol_dy3 = (x1*z2)/6 - (x2*z1)/6 - (x1*z4)/6 + (x4*z1)/6 + (x2*z4)/6 - (x4*z2)/6;
dvoly3_dp_i = vol_dy3 * dy3_pi;
vol_dz3 = (x2*y1)/6 - (x1*y2)/6 + (x1*y4)/6 - (x4*y1)/6 - (x2*y4)/6 + (x4*y2)/6;
dvolz3_dp_i = vol_dz3 * dz3_pi;

vol_dx4 = (y1*z2)/6 - (y2*z1)/6 - (y1*z3)/6 + (y3*z1)/6 + (y2*z3)/6 - (y3*z2)/6;
dvolx4_dp_i = vol_dx4 * dx4_pi;
vol_dy4 = (x2*z1)/6 - (x1*z2)/6 + (x1*z3)/6 - (x3*z1)/6 - (x2*z3)/6 + (x3*z2)/6;
dvoly4_dp_i = vol_dy4 * dy4_pi;
vol_dz4 = (x1*y2)/6 - (x2*y1)/6 - (x1*y3)/6 + (x3*y1)/6 + (x2*y3)/6 - (x3*y2)/6;
dvolz4_dp_i = vol_dz4 * dz4_pi;

dvol_dPfree = dvolx1_dp_i + dvoly1_dp_i + dvolz1_dp_i + ...
    dvolx2_dp_i + dvoly2_dp_i + dvolz2_dp_i + ...
    dvolx3_dp_i + dvoly3_dp_i + dvolz3_dp_i + ...
    dvolx4_dp_i + dvoly4_dp_i + dvolz4_dp_i;