function derKc_Fdes = compute_derKc_Fdes(invABwXi,G_inv_A,G_inv_ABw,Xi,derG_Fdes,derA_Fdes,derBw_Fdes,derXi_Fdes)
    derKc_Fdes = derG_Fdes * invABwXi - G_inv_A * derA_Fdes * invABwXi + G_inv_A * derBw_Fdes * Xi + G_inv_ABw * derXi_Fdes;
end