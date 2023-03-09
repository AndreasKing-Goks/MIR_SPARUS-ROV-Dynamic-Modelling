function squareprism_inertia = f(mass,spanx,spany,spanz)
    Ix = mass/12*(spany^2+spanz^2);
    Iy = mass/12*(spanx^2+spanz^2);
    Iz = mass/12*(spanx^2+spany^2);
    squareprism_inertia = [Ix 0 0;
                        0 Iy 0;
                        0 0 Iz];
end
