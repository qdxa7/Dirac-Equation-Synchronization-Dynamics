load('../SBM_M4_N200')

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

%%%%%%%%%%%%%%%% PLOTTING CODE STARTS BELOW %%%%%%%%%%%%%
figure
colors = lines(3);

%%% === TILE 1: Spectrum plot with inset and arrows ===
ax1 = gca;
hold(ax1, 'on');
plot(ax1, dnneglam, rho_nneg, 'o', 'Color', 'k');
set(ax1, 'FontSize', 16, 'TickLabelInterpreter', 'latex');
xlabel(ax1, '$E$', 'Interpreter','latex','FontSize',25,'FontWeight','bold');
ylabel(ax1, '$\rho_c(E)$', 'Interpreter','latex','FontSize',20,'FontWeight','bold');
ylim(ax1, [0 0.17]);
rectangle(ax1, 'Position', [1.28, 0.105, 0.4, 0.015], ...
         'EdgeColor', 'r', 'LineWidth', 1.5);
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