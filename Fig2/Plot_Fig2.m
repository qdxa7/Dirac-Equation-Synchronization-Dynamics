figure('Position', [70, 70, 1600, 900]) 
t = tiledlayout(2,4, 'Padding','compact', 'TileSpacing','compact');

%%Panel (a)
load('DESD_m1_E1_fig2.mat')

ax1 = nexttile(t,1);
p=plot(ax1,G);
p.NodeLabel = {};
p.NodeCData = estate(1:N,sel);
p.EdgeColor = 'k';
p.MarkerSize = 7;
p.LineWidth = 1;
p.EdgeAlpha = 1;
ax1.CLim = [min(estate(1:N,sel)), max(estate(1:N,sel))];
colormap(ax1,"turbo")
cb1 = colorbar(ax1,'westoutside');
cb1.Label.String = '$\Psi^{({\bar E})}_{i}$'; 
cb1.Label.Interpreter = 'latex';
cb1.FontSize = 11;
cb1.Label.FontSize = 25;
text(-0.50, 1, '(a)', 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 20,'Interpreter', 'latex')
axis(ax1,'off');

%%Panel (b)
ax2 = nexttile(t,2);
p=plot(ax2,G);
p.NodeLabel = {};
p.NodeColor = 'k';
p.EdgeCData = estate(N+1:end,sel);
p.MarkerSize = 1;
p.LineWidth = 3;
p.EdgeAlpha = 1;
ax2.CLim = [min(estate(N+1:end,sel)), max(estate(N+1:end,sel))];
colormap(ax2,"turbo")
cb2 = colorbar(ax2,'westoutside');
cb2.Ticks = linspace(-0.05, 0.05, 5); 
cb2.Label.String = '$\Psi^{({\bar E})}_{\ell}$'; 
cb2.Label.Interpreter = 'latex';
cb2.FontSize = 11;
cb2.Label.FontSize = 25;
text(-0.50, 1, '(b)', 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 20,'Interpreter', 'latex')
axis(ax2,'off');

%%Panel (c)
ax3 = nexttile(t,3);
hold(ax3, 'on');
scatter(Psi_dot_at_Tmax(1:N), estate(1:N,sel), 40, 'filled') 
set(gca, 'FontSize', 16, 'TickLabelInterpreter', 'latex');
xlabel('$\dot{\theta}_{i}$', 'Interpreter', 'latex','FontSize', 25)
ylabel('$\Psi^{({\bar E})}_{i}$', 'Interpreter', 'latex','FontSize', 25)
text(-0.35, 1, '(c)', 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 20,'Interpreter', 'latex')
box on

%%Panel (d)
ax4 = nexttile(t,4);
hold(ax4, 'on');
scatter(Psi_dot_at_Tmax(N+1:end), estate(N+1:end,sel), 40, 'filled')
set(gca, 'FontSize', 16, 'TickLabelInterpreter', 'latex');
xlabel('$\dot{\phi}_{\ell}$', 'Interpreter', 'latex','FontSize', 25)
ylabel('$\Psi^{({\bar E})}_{\ell}$', 'Interpreter', 'latex','FontSize', 25)
text(-0.40, 1, '(d)', 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 20,'Interpreter', 'latex')
box on

%%Panel (e)
load('DESD_m1_E2_fig2.mat')

ax5 = nexttile(t,5);
p=plot(ax5,G);
p.NodeLabel = {};
p.NodeCData = estate(1:N,sel);
p.EdgeColor = 'k';
p.MarkerSize = 7;
p.LineWidth = 1;
p.EdgeAlpha = 1;
ax5.CLim = [min(estate(1:N,sel)), max(estate(1:N,sel))];
colormap(ax5,"turbo")
cb5 = colorbar(ax5,'westoutside');
cb5.Label.String = '$\Psi^{({\bar E})}_{i}$';
cb5.Label.Interpreter = 'latex';
cb5.FontSize = 11;
cb5.Label.FontSize = 25;
text(-0.50, 1,'(e)', 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 20,'Interpreter', 'latex')
axis(ax5,'off');

%%Panel (f)
ax6 = nexttile(t,6);
p=plot(ax6,G);
p.NodeLabel = {};
p.NodeColor = 'k';
p.EdgeCData = estate(N+1:end,sel);
p.MarkerSize = 1;
p.LineWidth = 3;
p.EdgeAlpha = 1;
ax6.CLim = [min(estate(N+1:end,sel)), max(estate(N+1:end,sel))];
colormap(ax6,"turbo")
cb6 = colorbar(ax6,'westoutside');
cb6.Label.String = '$\Psi^{({\bar E})}_{\ell}$';
cb6.Label.Interpreter = 'latex';
cb6.FontSize = 11;
cb6.Label.FontSize = 25;
text(-0.50, 1, '(f)', 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 20,'Interpreter', 'latex')
axis(ax6,'off');

%%Panel (g)
ax7 = nexttile(t,7);
hold(ax7, 'on');
scatter(Psi_dot_at_Tmax(1:N), estate(1:N,sel), 40, 'filled')
set(gca, 'FontSize', 16, 'TickLabelInterpreter', 'latex');
xlabel('$\dot{\theta}_{i}$', 'Interpreter', 'latex','FontSize', 25)
lims = xlim(ax7);
ax7.XLim = [lims(1) 0.2]; 
xticks([-0.2 0 0.2])
ylabel('$\Psi^{({\bar E})}_{i}$', 'Interpreter', 'latex','FontSize', 25)
text(-0.35, 1, '(g)', 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 20,'Interpreter', 'latex')
box on

%%Panel (h)
ax8 = nexttile(t,8);
hold(ax8, 'on');
scatter(Psi_dot_at_Tmax(N+1:end), estate(N+1:end,sel), 40, 'filled')
set(gca, 'FontSize', 16, 'TickLabelInterpreter', 'latex');
xlabel('$\dot{\phi}_{\ell}$', 'Interpreter', 'latex','FontSize', 25)
ylabel('$\Psi^{({\bar E})}_{\ell}$', 'Interpreter', 'latex','FontSize', 25)
xticks([-0.05 0 0.05])
ax8.YTickMode = 'manual';
yticks([-0.05 0 0.05])
text(-0.40, 1, '(h)', 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 20,'Interpreter', 'latex')
box on