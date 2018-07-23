clear all;
close all;

%--------------------------------------------------------------------------
% Parameters

Parameters;

%--------------------------------------------------------------------------

[ k_rw_i, k_rn_i ] = relativePermeabilities( S, K_rw_i, K_rn_i, S_wr_i, S_nr_i, alpha_i, beta_i );
[ k_rw_d, k_rn_d ] = relativePermeabilities( S, K_rw_d, K_rn_d, S_wr_d, S_nr_d, alpha_d, beta_d );

f_i = f( S, K_rw_i, K_rn_i, S_wr_i, S_nr_i, alpha_i, beta_i, mu_w, mu_n ); 
f_d = f( S, K_rw_d, K_rn_d, S_wr_d, S_nr_d, alpha_d, beta_d, mu_w, mu_n ); 
                           
h_i = k_rn_i.*f_i;
h_d = k_rn_d.*f_d;

%--------------------------------------------------------------------------

figure(1)

subplot(1,3,1)
plot(S,k_rw_i,'r','LineWidth',3);
hold  on;
plot(S,k_rw_d,'b','LineWidth',3);
hold  off;
xlim([0.0 1.0-S_wr_d]);
xlabel('$S$','FontSize',20,'Interpreter','latex');
ylabel('$k_{rw}$','FontSize',20,'Interpreter','latex');
legend('k_{rw_i}','k_{rw_d}','Location','NorthWest');
grid  on;
set(gca,'FontSize',18);

subplot(1,3,2)
plot(S,k_rn_i,'r','LineWidth',3);
hold  on;
plot(S,k_rn_d,'b','LineWidth',3);
hold  off;
xlim([S_nr_d 1.0]);
xlabel('$S$','FontSize',20,'Interpreter','latex');
ylabel('$k_{rn}$','FontSize',20,'Interpreter','latex');
legend('k_{rn_i}','k_{rn_d}','Location','NorthEast');
grid  on;
set(gca,'FontSize',18);

subplot(1,3,3)
plot(S,h_i,'r','LineWidth',3);
hold  on;
plot(S,h_d,'b','LineWidth',3);
hold  off;
xlim([S_nr_d 1.0-S_wr_d]);
xlabel('$S$','FontSize',20,'Interpreter','latex');
ylabel('$h$','FontSize',20,'Interpreter','latex');
legend('h_i','h_d','Location','NorthWest');
grid  on;
set(gca,'FontSize',18);

%--------------------------------------------------------------------------

S_pc_i = (S_wr_i+0.0001):0.00001:(1.0-S_nr_i);
pc_i   = Pc(S_pc_i,Pb_i,S_wr_i,S_nr_i,gamma_i);

S_pc_d = (S_wr_d+0.0001):0.00001:(1.0-S_nr_d);
pc_d   = Pc(S_pc_d,Pb_d,S_wr_d,S_nr_d,gamma_d);

[ pc_plus, pc_minus ] = Pc_pm( S_pc_i, Pb_i, Pb_d, S_wr_i, S_wr_d,...
                               S_nr_i, S_nr_d, gamma_i, gamma_d );

%--------------------------------------------------------------------------

xi = (-1+0.001):0.000001:(1-0.001);

psi_eps = Psi( xi, epsilon );

figure(2)
plot(xi,psi_eps,'r','LineWidth',2);
hold on
line([1;1],[25;-25],'Color','k','LineStyle','--');
line([-1;-1],[25;-25],'Color','k','LineStyle','--');
hold off
xlabel('$\xi$','FontSize',20,'Interpreter','latex');
ylabel('$\Psi_{\epsilon}$','FontSize',20,'Interpreter','latex');
xlim([-1.05 1.05]);
ylim([-25 25]);
grid  on;
set(gca,'FontSize',18);

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
                             
                             
 G_plus_x = [];
 G_plus_y = [];
 
 G_minus_x = [];
 G_minus_y = [];
 
for S_g= 0.001:0.0025:0.99
    
  p_d = Pc(S_g,Pb_d,S_wr_d,S_nr_d,gamma_d);  
  p_i = Pc(S_g,Pb_i,S_wr_i,S_nr_i,gamma_i);
    
  p=p_i:1.0:p_d;
       
  G = G_function( S_g, p/P_ref, c, A, N_c,...
                  K_rw_i, K_rw_d, K_rn_i, K_rn_d,...
                  Pb_i/P_ref, Pb_d/P_ref, S_wr_i, S_wr_d, S_nr_i, S_nr_d,...
                  alpha_i, alpha_d, beta_i, beta_d, gamma_i, gamma_d,...
                  mu_w, mu_n );
              
  G_plus_x = [ G_plus_x;S_g*ones( length( p( G>0 ) ),1 ) ];
  G_plus_y = [ G_plus_y;p( G>0 )' ];

  G_minus_x = [ G_minus_x;S_g*ones( length( p( G<0 ) ),1 ) ];
  G_minus_y = [ G_minus_y;p( G<0 )' ];
  
end

[t,y] = ode23s(@(t,y) Rhs_ODE( t, y, c, A, N_c,...
                               K_rw_i, K_rw_d, K_rn_i, K_rn_d,...
                               Pb_i/P_ref, Pb_d/P_ref, S_wr_i, S_wr_d, S_nr_i, S_nr_d,...
                               alpha_i, alpha_d, beta_i, beta_d, gamma_i, gamma_d,...
                               mu_w, mu_n, epsilon ), [0 xi_end],...
                               [S_r+0.001;pc_plus_r/P_ref] );

                          
figure(3)
subplot(1,2,1)
plot(S_pc_i,pc_i,'r','LineWidth',2);
hold on
plot(S_pc_d,pc_d,'b','LineWidth',2);
plot(S_pc_i,pc_plus,'k--','LineWidth',2);
plot(S_l,pc_plus_l,'o','MarkerSize',10,'MarkerFaceColor','b');
plot(S_r,pc_plus_r,'o','MarkerSize',10,'MarkerFaceColor','r');
plot(G_plus_x,G_plus_y,'.','MarkerSize',1,'MarkerFaceColor','g');
plot(G_minus_x,G_minus_y,'.','MarkerSize',1,'MarkerFaceColor','k');
plot(y(:,1),y(:,2)*P_ref,'g--','LineWidth',2);
hold off
ylim([0 940]);
xlabel('$S$','FontSize',20,'Interpreter','latex');
ylabel('$p_c$','FontSize',20,'Interpreter','latex');
legend('pc_i','pc_d','pc_{plus}','E_l','E_r','Location','NorthEast');
grid  on;
set(gca,'FontSize',16);

subplot(1,2,2)
plot(S,f_i,'r','LineWidth',2);
hold  on;
plot(S,f_d,'b','LineWidth',2);
plot(S,0.5*(f_i+f_d),'k--','LineWidth',2);
plot(S,c*(S-S_r)+f_r,'g','LineWidth',2);
plot(S_l,f_l,'o','MarkerSize',10,'MarkerFaceColor','b');
plot(S_r,f_r,'o','MarkerSize',10,'MarkerFaceColor','r');
hold  off;
xlim([0.0 1.0]);
xlabel('$S$','FontSize',20,'Interpreter','latex');
ylabel('$f$','FontSize',20,'Interpreter','latex');
legend('f_i','f_d','f_p','c*S-A','E_l','E_r','Location','NorthWest');
grid  on;
set(gca,'FontSize',18);

