function eig_maxZ = matrixEigenValue(Z1,Z2,Z3,Z4,Z5,Z6,Z7,Z8,Z9)
Z = [Z1 Z2 Z3; Z4 Z5 Z6; Z7 Z8 Z9];
eig_maxZ = max(eig(Z));
end