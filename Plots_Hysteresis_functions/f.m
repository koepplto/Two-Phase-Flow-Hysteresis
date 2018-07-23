function [ f_value ] = f( S, K_rw, K_rn, S_wr, S_nr, alpha, beta, mu_w, mu_n )

[ k_rw, k_rn ] = relativePermeabilities( S, K_rw, K_rn, S_wr, S_nr, alpha, beta );

f_value = k_rw./( k_rw + (mu_w/mu_n)*k_rn );

end

