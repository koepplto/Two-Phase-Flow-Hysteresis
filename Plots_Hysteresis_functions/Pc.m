function [ p_c ] = Pc( S, Pb, S_wr, S_nr, alpha )

S_e = (S-S_wr)/(1-S_wr-S_nr);

p_c = Pb*(S_e.^(-1/alpha)-1).^(1-alpha);

end

