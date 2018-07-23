function [ h_h ] = h_hysteresis( S, p, K_rw_i, K_rw_d, K_rn_i, K_rn_d,...
                                 Pb_i, Pb_d, S_wr_i, S_wr_d, S_nr_i, S_nr_d,...
                                 alpha_i, alpha_d, beta_i, beta_d, gamma_i, gamma_d,...
                                 mu_w, mu_n )

 [ pc_p, pc_m ] = Pc_pm( S, Pb_i, Pb_d, S_wr_i, S_wr_d,...
                         S_nr_i, S_nr_d, gamma_i, gamma_d );

  h_i = h( S, K_rw_i, K_rn_i, S_wr_i, S_nr_i, alpha_i, beta_i, mu_w, mu_n );
  
  h_d = h( S, K_rw_d, K_rn_d, S_wr_d, S_nr_d, alpha_d, beta_d, mu_w, mu_n );
  
  h_p = 0.5*(h_i+h_d);
  
  h_m = 0.5*(h_d-h_i);
  
  h_h = h_p - h_m.*( (pc_p-p)./pc_m );
  
end

