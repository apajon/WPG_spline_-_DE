function derXi_Zdes = compute_derXi_Zdes(dOl_dZdes,index_Zdes,ind,contAngle)
    if ind==contAngle
        if index_Zdes==1
            derXi_Zdes = [0, 0, 0,      0,                  dOl_dZdes(3),           -dOl_dZdes(2);
                          0, 0, 0,      -dOl_dZdes(3),          0,                  dOl_dZdes(1)-1;
                          0, 0, 0,      dOl_dZdes(2),  -dOl_dZdes(1)+1, 0;
                          0, 0, 0, 0, 0, 0;
                          0, 0, 0, 0, 0, 0;
                          0, 0, 0, 0, 0, 0];
        else
            derXi_Zdes = [0, 0, 0,      0,                  dOl_dZdes(3),           -dOl_dZdes(2)+1;
                          0, 0, 0,      -dOl_dZdes(3),          0,                  dOl_dZdes(1);
                          0, 0, 0,      dOl_dZdes(2)-1,  -dOl_dZdes(1), 0;
                          0, 0, 0, 0, 0, 0;
                          0, 0, 0, 0, 0, 0;
                          0, 0, 0, 0, 0, 0];        
        end
    else
        derXi_Zdes = [0, 0, 0,      0,                  dOl_dZdes(3),           -dOl_dZdes(2);
                      0, 0, 0,      -dOl_dZdes(3),          0,                  dOl_dZdes(1);
                      0, 0, 0,      dOl_dZdes(2),  -dOl_dZdes(1), 0;
                      0, 0, 0, 0, 0, 0;
                      0, 0, 0, 0, 0, 0;
                      0, 0, 0, 0, 0, 0];            
    end
end