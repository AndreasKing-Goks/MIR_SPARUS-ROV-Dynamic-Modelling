function cylinder_inertia = f(mass,radius,height)
    Ix = (mass*radius^2)/2;
    Iy = mass*((radius^2)/4+(height^2)/12);
    Iz = Iy;
    cylinder_inertia = [Ix 0 0;
                        0 Iy 0;
                        0 0 Iz];
end