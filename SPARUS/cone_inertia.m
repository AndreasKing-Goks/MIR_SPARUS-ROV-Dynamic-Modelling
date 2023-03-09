function cone_inertia = f(mass,radius,height)
    Iz = 3/5*(mass*height^2)+3/20*(mass*radius^2);
    Ix = 3/10*mass*radius^2;
    Iy = Iz;
    cone_inertia = [Ix 0 0;
                    0 Iy 0;
                    0 0 Iz];
end