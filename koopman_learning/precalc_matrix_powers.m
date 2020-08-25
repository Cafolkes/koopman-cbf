function K_pows = precalc_matrix_powers(N_max,K)
    temp = K;
    for i=1:N_max
        K_pows{i}=temp;
        temp = temp*K;
    end