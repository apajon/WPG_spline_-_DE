function derKc_Zdes = compute_derKc_Zdes(invABwXi,G_inv_A,G_inv_ABw,Xi,derG_Zdes,derA_Zdes,derBw_Zdes,derXi_Zdes)
    derKc_Zdes = derG_Zdes * invABwXi - G_inv_A * derA_Zdes * invABwXi + G_inv_A * derBw_Zdes * Xi + G_inv_ABw * derXi_Zdes;
end