figure('Position', [70, 70, 1400, 900])
t = tiledlayout(2,3, 'Padding','compact', 'TileSpacing','compact');

cmap = lines(6);
colormap(cmap); 

%%Panel (a)
ax1 = nexttile(t,1);
hold(ax1, 'on');

load('DESD_m1_E1_fig9.mat')

theta_dot = Psi_dot_at_Tmax(1:N);
phi_dot = Psi_dot_at_Tmax(N+1:L);

gn = N/V; %number of groups
gl = repelem(1:gn,V);
gl = gl'; %group label

sel_estate = estate(1:N,sel);

p = plot(G);
p.NodeCData = gl;
p.EdgeColor = 'k';
p.MarkerSize = 10;
p.LineWidth = 1;
p.EdgeAlpha = 1;
p.NodeLabel = {};
text(0, 1, '(a)', 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 20,'Interpreter', 'latex')
axis off


%%Panel (b)
ax2 = nexttile(t,2);
hold(ax2, 'on');

scatter(Omega(1:N), sel_estate, 40, gl, 'filled')
set(gca, 'FontSize', 16, 'TickLabelInterpreter', 'latex');
xlabel('$\omega_{i}$', 'Interpreter', 'latex','FontSize', 25)
ylabel('$\Psi^{({\bar E})}_{i}$', 'Interpreter', 'latex','FontSize', 25)
xlim([-4 4])
text(-0.25, 1, '(b)', 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 20,'Interpreter', 'latex')
box on


%%Panel (c)
ax3 = nexttile(t,3);
hold(ax3, 'on');

scatter(Omega(1:N), theta_dot, 40, gl, 'filled') 
set(gca, 'FontSize', 16, 'TickLabelInterpreter', 'latex');
xlabel('$\omega_{i}$', 'Interpreter', 'latex','FontSize', 25)
ylabel('$\dot{\theta}_{i}$', 'Interpreter', 'latex','FontSize', 25)
xlim([-4 4])
text(-0.2, 1, '(c)', 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 20,'Interpreter', 'latex')
box on


%%Panel (d)
ax4 = nexttile(t,4);
hold(ax4, 'on');

gn = N/V; %number of groups
gl = repelem(1:gn,V);
gl = gl'; %group label

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

sel_estate = estate(N+1:N+L,sel);
cscale = turbo(100);

p = plot(G);
p.NodeColor = 'k';
p.EdgeCData = edge_colors;
p.EdgeColor = 'flat';
p.MarkerSize = 8;
p.LineWidth = 2;
p.EdgeAlpha = 1;
p.NodeLabel = {};
text(0, 1, '(d)', 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 20,'Interpreter', 'latex')
axis off


%%Panel (e)
ax5 = nexttile(t,5);
hold(ax5, 'on');

scatter(Omega(inter_idx), sel_estate(inter_idx), 40, new_labels, 'filled')
set(gca, 'FontSize', 16, 'TickLabelInterpreter', 'latex');
xlabel('$\hat{\omega}_{\ell}$', 'Interpreter', 'latex','FontSize', 25)
ylabel('$\Psi^{(\bar{E})}_{\ell}$', 'Interpreter', 'latex','FontSize', 25)
xlim([-4 4])
text(-0.25, 1, '(e)', 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 20,'Interpreter', 'latex')
box on


%%Panel (f)
ax6 = nexttile(t,6);
hold(ax6, 'on');

scatter(Omega(inter_idx), phi_dot(inter_idx), 40, new_labels, 'filled')
set(gca, 'FontSize', 16, 'TickLabelInterpreter', 'latex');
xlabel('$\hat{\omega}_{\ell}$', 'Interpreter', 'latex','FontSize', 25)
ylabel('$\dot{\phi}_{\ell}$', 'Interpreter', 'latex','FontSize', 25)
xlim([-4 4])
text(-0.2, 1, '(f)', 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 20,'Interpreter', 'latex')
box on