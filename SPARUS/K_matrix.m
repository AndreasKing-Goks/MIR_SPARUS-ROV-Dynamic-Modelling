function K_matrix = f(L, Sx, Sy, Sz, Dy, Dz, Cd11, Cd22, Cd33)
    rho = 1000;

    k11 = rho*Cd11*Sx/2;
    k22 = rho*Cd22*Sy/2;
    k33 = rho*Cd22*Sz/2;
    k44 = 0;
    k55 = rho*Cd33*Dz*(L^4)/64;
    k66 = rho*Cd33*Dy*(L^4)/64;

    K_matrix = [k11 0 0 0 0 0 ;
                   0 k22 0 0 0 0 ;
                   0 0 k33 0 0 0 ;
                   0 0 0 k44 0 0 ;
                   0 0 0 0 k55 0 ;
                   0 0 0 0 0 k66];
end