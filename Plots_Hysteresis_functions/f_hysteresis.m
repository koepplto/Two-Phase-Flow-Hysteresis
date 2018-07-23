function [ f_h ] = f_hysteresis( S, p, K_rw_i, K_rw_d, K_rn_i, K_rn_d,...
                                 Pb_i, Pb_d, S_wr_i, S_wr_d, S_nr_i, S_nr_d,...
                                 alpha_i, alpha_d, beta_i, beta_d, gamma_i, gamma_d,...
                                 mu_w, mu_n )

 [ pc_p, pc_m ] = Pc_pm( S, Pb_i, Pb_d, S_wr_i, S_wr_d,...
                         S_nr_i, S_nr_d, gamma_i, gamma_d );

  f_i = f( S, K_rw_i, K_rn_i, S_wr_i, S_nr_i, alpha_i, beta_i, mu_w, mu_n );
  
  f_d = f( S, K_rw_d, K_rn_d, S_wr_d, S_nr_d, alpha_d, beta_d, mu_w, mu_n );
  
  f_p = 0.5*(f_i+f_d);
  
  f_m = 0.5*(f_d-f_i);
  
  f_h = f_p - f_m.*((pc_p-p)./pc_m);
  
end

