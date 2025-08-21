% Sample input (replace with your actual data)
% dnneglam, rho_nneg, neglam, rho_neg, dE, pidx, lowlam, hilam, N must be defined

tic

load("../RG_k12_N1000");

%%Generate Dirac operator
D=sparse(N+L,N+L);
D(1:N, N+1:N+L)=B1;
D(N+1:N+L, 1:N)=transpose(B1);
gamma=sparse(N+L,N+L);
gamma(1:N, 1:N)=speye(N);
gamma(N+1:N+L, N+1:N+L)=-speye(L);
m = 1; %specify the mass if interested in the Dirac operator.
Dirac = D + (m*gamma);

%%Caluclate spectra and sort
[eve,eva] = eigs(Dirac,length(Dirac));
lam = diag(eva);
[slam,idx2] = sort(lam,'ascend'); %sort eigenvalues from lowest to highest
seve = eve(:,idx2); %sort eigenvectors compatibly with eigenvalues

%%Split eigenvalues into negative and non-negative values
[~,~,neglam]=find(slam.*(slam<0));
[idx3,~]=find(slam>=0);
nneglam = slam(idx3);

dnneglam = flip(nneglam);
rho_nneg = [1:numel(nneglam)]/numel(lam);
rho_neg = [1:numel(neglam)]/numel(lam);

pidx = 500;
E = nneglam(pidx);

%%%%%%%%%%%%%%%% PLOTTING CODE STARTS BELOW %%%%%%%%%%%%%

figure('Position', [100, 100, 1300, 1000])  % wide figure
t = tiledlayout(2,2, 'Padding','compact', 'TileSpacing','compact');

%%% === TILE 1: Spectrum plot with inset and arrows ===
ax1 = nexttile(t,1);
hold(ax1, 'on');
% axis(ax1, 'square')  % force square plot
plot(ax1, dnneglam, rho_nneg, 'o', 'Color', 'k');
set(ax1, 'FontSize', 16, 'TickLabelInterpreter', 'latex');
xlabel(ax1, '$E$', 'Interpreter','latex','FontSize',25,'FontWeight','bold');
ylabel(ax1, '$\rho_c(E)$', 'Interpreter','latex','FontSize',20,'FontWeight','bold');
ylim(ax1, [0 0.17]);
text(-0.13, 1, '(a)', 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 20,'Interpreter', 'latex')
box(ax1, 'on');

%%% === INSET AXES ===
main_pos = ax1.Position;
inset_pos = [main_pos(1)+0.58*main_pos(3), ... % horizontal position
    main_pos(2)+0.59*main_pos(4), ... % vertical position
    0.38*main_pos(3), ...             % width
    0.38*main_pos(4)];                % height
inset_ax = axes('Position', inset_pos);
box(inset_ax, 'on'); hold(inset_ax, 'on');
plot(inset_ax, dnneglam, rho_nneg, '.', 'Color', 'k');
plot(inset_ax, neglam, rho_neg, '.', 'Color', 'k');
rectangle(inset_ax, 'Position', [0.1, 0, 6, 0.3], ...
    'EdgeColor', 'r', 'LineWidth', 1.5);
xlabel(inset_ax, '$E$', 'Interpreter', 'latex', 'FontSize', 16);
ylabel(inset_ax, '$\rho_c(E)$', 'Interpreter', 'latex', 'FontSize', 16);
set(inset_ax, 'FontSize', 12, 'TickLabelInterpreter', 'latex');
%%% === ARROWS (annotations must be in figure units) ===
ax_pos = get(ax1, 'Position');
colors = get(gca, 'ColorOrder');
% Arrow 1: Convert data point to normalized figure units
px1 = ((nneglam(2) - ax1.XLim(1)) / diff(ax1.XLim) * ax_pos(3) + ax_pos(1))+0.003;
py1 = ((rho_nneg(N - (2 - 1)) - ax1.YLim(1)) / diff(ax1.YLim) * ax_pos(4) + ax_pos(2))+0.01;
% Define arrow start (a bit away from point)
px_start1 = px1 + 0.0;  % shift right
py_start1 = py1 + 0.03;  % shift up
% Draw arrow
annotation('textarrow', [px_start1, px1], [py_start1, py1], 'Color', colors(1,:)); %'String', sprintf('E = %0.4f',E));
% Arrow 2: Convert data point to normalized figure units
px2 = (nneglam(500) - ax1.XLim(1)) / diff(ax1.XLim) * ax_pos(3) + ax_pos(1);
py2 = ((rho_nneg(N - (500 - 1)) - ax1.YLim(1)) / diff(ax1.YLim) * ax_pos(4) + ax_pos(2))+0.035;
% Define arrow start (a bit away from point)
px_start2 = px2 + 0.0;  % shift right
py_start2 = py2 + 0.03;  % shift up
% Draw arrow
annotation('textarrow', [px_start2, px2], [py_start2, py2], 'Color', colors(2,:)); %'String', sprintf('E = %0.4f',E));


%%% === TILE 2: Log-log slope plot ===
ax2 = nexttile(t,2);
hold(ax2, 'on');
% axis(ax2, 'square')  % force square plot

lowlam = 50;
hilam = 50;
dE = abs(nneglam - nneglam(pidx));
x1 = log(dE(pidx-1:-1:lowlam));
x2 = log(dE(pidx+1:pidx+hilam));
l_rho = (1:numel(x1))/numel(x1);
g_rho = (1:numel(x2))/numel(x2);

plot(ax2, x1, log(l_rho), 'o', 'Color', colors(2,:), 'DisplayName','$\bar{E}-E < 0$');
plot(ax2, x2, log(g_rho), '*', 'Color', colors(2,:), 'DisplayName','$\bar{E}-E > 0$');

% Reference slope = 1 lines
y_ref1 = x1 - x1(end) + log(l_rho(end))+0.25;
y_ref2 = x2 - x2(end) + log(g_rho(end))-0.1;
plot(ax2, x1, y_ref1, '--', 'LineWidth', 1, 'Color', 'k', 'DisplayName', '$\delta=1$');
plot(ax2, x2, y_ref2, '--', 'LineWidth', 1, 'Color', 'k', 'HandleVisibility','off');

set(ax2, 'FontSize', 16, 'TickLabelInterpreter', 'latex');
xlabel(ax2, 'ln($|\bar{E}-E|$)', 'Interpreter','latex','FontSize',20,'FontWeight','bold');
ylabel(ax2, 'ln($\rho_c(E)$)', 'Interpreter','latex','FontSize',20,'FontWeight','bold');
legend(ax2,'Box','off','Interpreter','latex','FontSize',20);
text(-0.13, 1, '(b)', 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 20,'Interpreter', 'latex')
box(ax2, 'on');

%%% === TILE 3
Elist = [2 4];
ax3 = nexttile(t,3);
hold(ax3, 'on');
for i = 1:2
    load(sprintf('DS_Dyn2s_m1_E%d_lam1_0_RG_k12_N1000_NPsi0',Elist(i)),'E','Sigma','sel2','seve','Psi_dot_fore','Psi_fore')

    pPsi_dot_fore = seve' * Psi_dot_fore;
    pPsi_dot_fore_full = pPsi_dot_fore;
    pPsi_dot_fore_sel = pPsi_dot_fore(sel2,:);
    pPsi_dot_fore(sel2,:) = [];
    V2 = mean(pPsi_dot_fore.^2);

    plot(ax3,Sigma,sqrt(V2),'linewidth', 2,'displayname', sprintf('$\\bar{E}$=%0.4f',E));
    legend('location',' best','Interpreter', 'latex','box','off','fontsize', 20,'fontweight','bold')
    set(gca, 'XScale', 'linear', 'YScale', 'log')

end
set(ax3, 'FontSize', 16, 'TickLabelInterpreter', 'latex');
xlabel('$\sigma$','Interpreter','latex','FontSize',25);
ylabel('$V$','Interpreter','latex','FontSize',20);
ylim([1e-6 1])
text(-0.13, 1, '(c)', 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 20,'Interpreter', 'latex')
box on

%%% === TILE 4
Elist = [2 4];
ax4 = nexttile(t,4);
hold(ax4, 'on');
for i = 1:2
    load(sprintf('DS_Dyn2t_m1_E%d_lam1_0_RG_k12_N1000',Elist(i)),'E','time','W2'); %,'sel2','seve','Psi_dot_fore','Psi_fore','time')

    plot(ax4,time,sqrt(W2),'linewidth', 2,'displayname', sprintf('$\\bar{E}$=%0.4f',E));
    legend('location',' best','Interpreter', 'latex','box','off','fontsize', 20,'fontweight','bold')
    set(gca, 'XScale', 'linear', 'YScale', 'log')

end
set(ax4, 'FontSize', 16, 'TickLabelInterpreter', 'latex');
xlabel('$t$','Interpreter','latex','FontSize',25);
ylabel('$W$','Interpreter','latex','FontSize',20);
text(-0.13, 1, '(d)', 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 20,'Interpreter', 'latex')
box on
