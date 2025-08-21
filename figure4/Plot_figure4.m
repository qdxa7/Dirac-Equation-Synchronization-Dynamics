%%% This plots Ralpha, Rbeta and F (the free energy) vs sigma

figure('Position', [100, 100, 500, 900])  % wide figure
t = tiledlayout(2,1, 'Padding','compact', 'TileSpacing','compact');

%%% === TILE 1
ax1 = nexttile(t,1);
hold(ax1, 'on');

load('DS_Dyn1s_m1_E2_lam1_0_RG_k12_N1000.mat','R_alpha_fore','R_beta_fore','Sigma')
plot(Sigma, R_alpha_fore,'linewidth', 2,'markersize', 5,'displayname', '$R_{\alpha}$'); hold on;
plot(Sigma, R_beta_fore,'linewidth', 2,'markersize', 5,'displayname', '$R_{\beta}$'); hold on;

set(gca, 'FontSize', 16, 'TickLabelInterpreter', 'latex');
xlabel('$\sigma$','Interpreter', 'latex','fontsize', 25, 'fontweight','bold')
legend('location','northeast','Interpreter', 'latex','box','off','fontsize', 16,'fontweight','normal')
text(-0.18, 1, '(a)', 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 20,'Interpreter', 'latex')
box on

%%% === TILE 2
ax2 = nexttile(t,2);
hold(ax2, 'on');

load('DS_Dyn1s_m1_E2_lam1_0_RG_k12_N1000.mat','Deqn','Psi_fore','Sigma')
modF = sum(cos(Deqn*Psi_fore));
normF = -modF/(size(Psi_fore,1));
F = Sigma.*modF;
plot(Sigma, normF,'linewidth', 2,'markersize', 5); hold on;

set(gca, 'FontSize', 16, 'TickLabelInterpreter', 'latex');
xlabel('$\sigma$','Interpreter', 'latex','fontsize', 25, 'fontweight','bold')
ylabel('$f$','Interpreter', 'latex','fontsize', 25, 'fontweight','normal')
text(-0.18, 1, '(b)', 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 20,'Interpreter', 'latex')
box on