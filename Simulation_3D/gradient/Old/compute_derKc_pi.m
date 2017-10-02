function derKc_p_i = compute_derKc_pi(invABwXi,G_inv_A,G_inv_ABw,Xi,derG_p_i,derA_p_i,derBw_p_i,derXi_p_i)
    derKc_p_i = derG_p_i * invABwXi - G_inv_A * derA_p_i * invABwXi + G_inv_A * derBw_p_i * Xi + G_inv_ABw * derXi_p_i;
end