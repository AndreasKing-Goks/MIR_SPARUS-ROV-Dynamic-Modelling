function halfsphere_inertia = f(mass,radius,height)
    Ix = 2/5*mass*radius^2;
    Iz = .259*mass*radius^2;
    Iy = Iz;
    halfsphere_inertia = [Ix 0 0;
                        0 Iy 0;
                        0 0 Iz];
end