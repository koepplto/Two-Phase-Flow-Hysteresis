function [ kr_w, kr_n ] = relativePermeabilities( S, K_rw, K_rn, S_wr, S_nr, alpha, beta )

S_e = (S-S_wr)/(1-S_wr-S_nr);

kr_w = K_rw*(S_e.^(1/2)).*( 1-( 1-S_e.^(1/alpha) ).^(alpha) ).^2;

kr_n = K_rn*(1-S_e).^(1/2).*( 1-S_e.^(1/beta) ).^(2*beta);

%kr_w = K_rw*S_e.^(2*alpha);

%kr_n = K_rn*(1-S_e.^(2*alpha));
end

