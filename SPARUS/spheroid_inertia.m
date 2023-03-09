function spheroid_inertia = f(mass,a,b)
    Ix = mass/5*2*b^2;
    Iy = mass/5*(a^2+b^2);
    Iz = Iy;
    spheroid_inertia = [Ix 0 0;
                        0 Iy 0;
                        0 0 Iz];
end
