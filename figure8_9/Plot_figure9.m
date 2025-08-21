%%% This plots the SBM network with theta and phi vs intrinsic frequency

%%% Tiled plots
figure('Position', [70, 70, 1600, 900])  % wide figure
t = tiledlayout(2,4, 'Padding','compact', 'TileSpacing','compact');


% === SET CUSTOM COLORMAP ===
cmap = lines(6);
colormap(cmap);  % Apply globally


%%% === TILE 1
ax1 = nexttile(t,1);
hold(ax1, 'on');
load('DS_Dyn1t_m1_E1_lam1_0_SBM_M4_N200_final1.mat')

% group label
gn = N/V; %number of groups
gl = repelem(1:gn,V);
gl = gl';

p = plot(G);
p.NodeCData = gl;
p.EdgeColor = 'k';
p.MarkerSize = 10;
p.LineWidth = 1;
p.EdgeAlpha = 1;
p.NodeLabel = {};
text(-0.12, 1, '(a)', 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 20,'Interpreter', 'latex')
axis off


%%% === TILE 2
ax2 = nexttile(t,2);
hold(ax2, 'on');

time_step = length(time);
theta_dot = Psi_dot_tot(1:N,time_step);

scatter(Omega(1:N), theta_dot, 40, gl, 'filled')  % 40 is point size
set(gca, 'FontSize', 16, 'TickLabelInterpreter', 'latex');
xlabel('$\omega_{i}$', 'Interpreter', 'latex','FontSize', 25)
ylabel('$\dot{\theta}_{i}$', 'Interpreter', 'latex','FontSize', 25)
xlim([-4 4])
text(-0.27, 1, '(b)', 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 20,'Interpreter', 'latex')
box on


%%% === TILE 3
ax3 = nexttile(t,3);
hold(ax3, 'on');
load('DS_Dyn1t_m1_E1_lam1_0_SBM_M4_N200_final2.mat','time','Omega','Psi_dot_tot')

time_step = length(time);
theta_dot = Psi_dot_tot(1:N,time_step);

scatter(Omega(1:N), theta_dot, 40, gl, 'filled')  % 40 is point size
set(gca, 'FontSize', 16, 'TickLabelInterpreter', 'latex');
xlabel('$\omega_{i}$', 'Interpreter', 'latex','FontSize', 25)
ylabel('$\dot{\theta}_{i}$', 'Interpreter', 'latex','FontSize', 25)
xlim([-4 4])
text(-0.27, 1, '(c)', 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 20,'Interpreter', 'latex')
box on


%%% === TILE 4
ax4 = nexttile(t,4);
hold(ax4, 'on');
load('DS_Dyn1t_m1_E2_lam1_0_SBM_M4_N200_final3.mat','time','Omega','Psi_dot_tot')

time_step = length(time);
theta_dot = Psi_dot_tot(1:N,time_step);

scatter(Omega(1:N), theta_dot, 40, gl, 'filled')  % 40 is point size
set(gca, 'FontSize', 16, 'TickLabelInterpreter', 'latex');
xlabel('$\omega_{i}$', 'Interpreter', 'latex','FontSize', 25)
ylabel('$\dot{\theta}_{i}$', 'Interpreter', 'latex','FontSize', 25)
xlim([-4 4])
text(-0.27, 1, '(d)', 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 20,'Interpreter', 'latex')
box on


%%% === TILE 5
ax5 = nexttile(t,5);
hold(ax5, 'on');
load('DS_Dyn1t_m1_E-1_lam1_0_SBM_M4_N200_final1.mat')

% group label
gn = N/V; %number of groups
gl = repelem(1:gn,V);
gl = gl';

% inter-community link indices
edge_list = G.Edges.EndNodes;
inter_idx = find(gl(edge_list(:,1)) ~= gl(edge_list(:,2)));
% create 'bundle' labels for links (equivalent to community labels for
% nodes)
% Get community pair for each edge
comm_pairs = [gl(edge_list(:,1)), gl(edge_list(:,2))];
% Sort each row so that (2,1) and (1,2) are treated as the same
comm_pairs_sorted = sort(comm_pairs, 2);  % ensures undirected symmetry
% Assign unique labels to each unique community pair
[~, ~, link_labels] = unique(comm_pairs_sorted, 'rows');
inter_link_labels = link_labels(inter_idx);
% Step 1: Get the unique values in inter_link_labels (in sorted order)
unique_vals = unique(inter_link_labels);  % returns [2 3 4 6 7 9]
% Step 2: Map each value to a new index 1 through 6 (for plot colouring purposes)
[~, new_labels] = ismember(inter_link_labels, unique_vals);
edge_colors = nan(L, 1);
edge_colors(inter_idx) = new_labels;  % assign only inter-community


p = plot(G);
p.NodeColor = 'k';
p.EdgeCData = edge_colors;
p.EdgeColor = 'flat';           % ensures coloring from EdgeCData is used
p.MarkerSize = 8;
p.LineWidth = 2;
p.EdgeAlpha = 1;
p.NodeLabel = {};
text(-0.12, 1,'(e)', 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 20,'Interpreter', 'latex')
axis off


%%% === TILE 6
ax6 = nexttile(t,6);
hold(ax6, 'on');

time_step = length(time);
phi_dot = Psi_dot_tot(N+1:L,time_step);

scatter(Omega(inter_idx), phi_dot(inter_idx), 40, new_labels, 'filled')  % 40 is point size
set(gca, 'FontSize', 16, 'TickLabelInterpreter', 'latex');
xlabel('$\hat{\omega}_{\ell}$', 'Interpreter', 'latex','FontSize', 25)
ylabel('$\dot{\phi}_{\ell}$', 'Interpreter', 'latex','FontSize', 25)
xlim([-4 4])
text(-0.27, 1, '(f)', 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 20,'Interpreter', 'latex')
box on


%%% === TILE 7
ax7 = nexttile(t,7);
hold(ax7, 'on');
load('DS_Dyn1t_m1_E-1_lam1_0_SBM_M4_N200_final2.mat','time','Omega','Psi_dot_tot')

time_step = length(time);
phi_dot = Psi_dot_tot(N+1:L,time_step);

scatter(Omega(inter_idx), phi_dot(inter_idx), 40, new_labels, 'filled')  % 40 is point size
set(gca, 'FontSize', 16, 'TickLabelInterpreter', 'latex');
xlabel('$\hat{\omega}_{\ell}$', 'Interpreter', 'latex','FontSize', 25)
ylabel('$\dot{\phi}_{\ell}$', 'Interpreter', 'latex','FontSize', 25)
xlim([-4 4])
text(-0.27, 1, '(g)', 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 20,'Interpreter', 'latex')
box on


%%% === TILE 8
ax8 = nexttile(t,8);
hold(ax8, 'on');
load('DS_Dyn1t_m1_E-2_lam1_0_SBM_M4_N200_final3.mat','time','Omega','Psi_dot_tot')

time_step = length(time);
phi_dot = Psi_dot_tot(N+1:L,time_step);

scatter(Omega(inter_idx), phi_dot(inter_idx), 40, new_labels, 'filled')  % 40 is point size
set(gca, 'FontSize', 16, 'TickLabelInterpreter', 'latex');
xlabel('$\hat{\omega}_{\ell}$', 'Interpreter', 'latex','FontSize', 25)
ylabel('$\dot{\phi}_{\ell}$', 'Interpreter', 'latex','FontSize', 25)
xlim([-4 4])
text(-0.27, 1, '(h)', 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 20,'Interpreter', 'latex')
box on