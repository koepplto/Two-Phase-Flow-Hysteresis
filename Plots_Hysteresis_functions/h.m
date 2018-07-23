function [ h_value ] = h( S, K_rw, K_rn, S_wr, S_nr, alpha, beta, mu_w, mu_n )

[ ~, k_rn ] = relativePermeabilities( S, K_rw, K_rn, S_wr, S_nr, alpha, beta );

f_value = f( S, K_rw, K_rn, S_wr, S_nr, alpha, beta, mu_w, mu_n );

h_value = k_rn.*f_value;

end

