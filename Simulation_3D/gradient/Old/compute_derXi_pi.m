function derXi_p_i = compute_derXi_pi(dOl_dp_i)
    derXi_p_i = [0, 0, 0,      0,                  dOl_dp_i(3),           -dOl_dp_i(2);
                0, 0, 0,      -dOl_dp_i(3),          0,                  dOl_dp_i(1);
                0, 0, 0,      dOl_dp_i(2),  -dOl_dp_i(1), 0;
                0, 0, 0, 0, 0, 0;
                0, 0, 0, 0, 0, 0;
                0, 0, 0, 0, 0, 0];   
end