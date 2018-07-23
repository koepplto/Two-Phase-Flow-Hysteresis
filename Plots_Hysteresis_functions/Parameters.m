
S = 0.0001:0.01:1.0-0.0001;

% Pb_i [pa]
Pb_i = 350.0;

% Pb_d [pa]
Pb_d = 700.0;

% S_wr_i [-]
S_wr_i = 0.0;

% S_wr_d [-]
S_wr_d = 0.0;

% S_nr_i [-]
S_nr_i = 0.0;

% S_wr_d [-]
S_nr_d = 0.0;

alpha_i = 0.6;
alpha_d = 0.85;

K_rw_i = 0.6;
K_rn_i = 1.0;

K_rw_d = 0.6;
K_rn_d = 0.6;

beta_i = 1.2;
beta_d = 1.1;

gamma_i = 0.9206;
gamma_d = 0.906;

% mu_w [Pas]
mu_w = 1.0e-3;

% mu_n [Pas]
mu_n = 3.0e-5;

S_l = 0.67;
S_r = 0.075;

% K [m^2]
K = 10.0^(-10);

% P_ref [Pa]
P_ref = 100.0;

% L_ref [m]
L_ref = 1.0;

% v_tot [m/s]
v_tot = 1.0e-5;

N_c = (K*P_ref)/(L_ref*v_tot*mu_n);
N_g = 0.0;

xi_end = 50.00;

epsilon = 0.01;
