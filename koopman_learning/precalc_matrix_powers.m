function [K_pows, CK_pows] = precalc_matrix_powers(N_max,K,C)
    temp = K;
    for i=1:N_max
        K_pows{i}=temp;
        CK_pows{i}=C*temp;
        temp = temp*K;
    end