function HT_matrix = f(p,q,r)
    S_matrix = [0 -r q;
                r 0 -p;
               -q p 0];

    HT_matrix = [eye(3) zeros(3);
                S_matrix eye(3)];
end