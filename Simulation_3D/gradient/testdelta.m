function d = testdelta(deltat,Fc_mc_splipt,Fc_mc_splipn)
fric = 0.8;
norm_deltat = norm(deltat);
d = norm_deltat * Fc_mc_splipt + fric * Fc_mc_splipn * deltat;
