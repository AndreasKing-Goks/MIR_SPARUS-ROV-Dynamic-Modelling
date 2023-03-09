function H_matrix = f(p,q,r)
    ST_matrix = [0 r -q;
                -r 0 p;
                q -p 0];

    H_matrix = [eye(3) ST_matrix;
                zeros(3) eye(3)];
end