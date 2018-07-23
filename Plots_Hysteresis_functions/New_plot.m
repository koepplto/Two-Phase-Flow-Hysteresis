clear all;
close all;

%--------------------------------------------------------------------------
% Parameters

Parameters;

%--------------------------------------------------------------------------

f_i_l = f( S_l, K_rw_i, K_rn_i, S_wr_i, S_nr_i, alpha_i, beta_i, mu_w, mu_n ); 
f_d_l = f( S_l, K_rw_d, K_rn_d, S_wr_d, S_nr_d, alpha_d, beta_d, mu_w, mu_n ); 

f_l = 0.5*(f_d_l+f_i_l);

f_i_r = f( S_r, K_rw_i, K_rn_i, S_wr_i, S_nr_i, alpha_i, beta_i, mu_w, mu_n ); 
f_d_r = f( S_r, K_rw_d, K_rn_d, S_wr_d, S_nr_d, alpha_d, beta_d, mu_w, mu_n ); 

f_r = 0.5*(f_d_r+f_i_r);

c = (f_r-f_l)/(S_r-S_l);

A = c*S_r-f_r;

[ pc_plus_l, pc_minus_l ] = Pc_pm( S_l, Pb_i, Pb_d, S_wr_i, S_wr_d,...
                                  S_nr_i, S_nr_d, gamma_i, gamma_d );

[ pc_plus_r, pc_minus_r ] = Pc_pm( S_r, Pb_i, Pb_d, S_wr_i, S_wr_d,...
                                 S_nr_i, S_nr_d, gamma_i, gamma_d );

%--------------------------------------------------------------------------

data = [];

epsilon=[0.75 0.5 0.1 0.02 0.005];

colors =['c' 'r--','b','k--','g'];

liste=cell(7,1);

for i=1:5

[t,y] = ode23s(@(t,y) Rhs_ODE( t, y, c, A, N_c,...
                               K_rw_i, K_rw_d, K_rn_i, K_rn_d,...
                               Pb_i/P_ref, Pb_d/P_ref, S_wr_i, S_wr_d, S_nr_i, S_nr_d,...
                               alpha_i, alpha_d, beta_i, beta_d, gamma_i, gamma_d,...
                               mu_w, mu_n, epsilon(i) ), [0 xi_end],...
                               [S_r+0.001;pc_plus_r/P_ref]);

liste{i} = ['\epsilon = ' num2str(epsilon(i))];

figure(1)

hold on
plot(t,y(:,1),colors(i),'LineWidth',2);
xlim([0 xi_end]);
ylim([0 1.0])
xlabel('$\xi$','FontSize',20,'Interpreter','latex');
ylabel('$S$','FontSize',20,'Interpreter','latex');
grid  on;
set(gca,'FontSize',18);

end

liste{6} = 'S_l';
liste{7} = 'S_r';

plot(0:0.01:xi_end,S_l*ones(length(0:0.01:xi_end),1),'k--','LineWidth',2);
plot(0:0.01:xi_end,S_r*ones(length(0:0.01:xi_end),1),'b--','LineWidth',2);
legend(liste,'Location','SouthEast');
hold off;
