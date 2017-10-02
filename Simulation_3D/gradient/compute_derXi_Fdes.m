function derXi_Fdes = compute_derXi_Fdes(dOl_dFdes)
    derXi_Fdes = [0, 0, 0,      0,                  dOl_dFdes(3),           -dOl_dFdes(2);
                0, 0, 0,      -dOl_dFdes(3),          0,                  dOl_dFdes(1);
                0, 0, 0,      dOl_dFdes(2),  -dOl_dFdes(1), 0;
                0, 0, 0, 0, 0, 0;
                0, 0, 0, 0, 0, 0;
                0, 0, 0, 0, 0, 0];   
end