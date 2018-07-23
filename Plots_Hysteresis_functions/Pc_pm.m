function [ pc_p, pc_m ] = Pc_pm( S, Pb_i, Pb_d, S_wr_i, S_wr_d,...
                                 S_nr_i, S_nr_d, gamma_i, gamma_d )

pc_d   = Pc(S,Pb_d,S_wr_d,S_nr_d,gamma_d);
pc_i   = Pc(S,Pb_i,S_wr_i,S_nr_i,gamma_i);

pc_p = 0.5*(pc_i+pc_d);
pc_m = 0.5*(pc_d-pc_i);

end

