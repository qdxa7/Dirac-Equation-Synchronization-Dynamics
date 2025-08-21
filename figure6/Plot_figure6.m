% Sample input (replace with your actual data)
% dnneglam, rho_nneg, neglam, rho_neg, dE, pidx, lowlam, hilam, N must be defined

figure('Position', [100, 100, 1300, 1300])  % wide figure
t = tiledlayout(2,2, 'Padding','compact', 'TileSpacing','compact');

%%% === TILE 1
ax1 = nexttile(t,1);
hold(ax1, 'on');
load('DS_Dyn1s_m1_E2_lam1_0_RG_k12_N1000.mat','E','sel2','slam','seve','Psi_dot_fore','Sigma')

pPsi_dot_fore = seve' * Psi_dot_fore;
pPsi_dot_fore_full = pPsi_dot_fore;
pPsi_dot_fore(sel2,:) = [];

rslam = slam;
nidx = find(slam<0); %if you want to remove E'<0 from plot
pPsi_dot_fore_alpha = pPsi_dot_fore_full;
pPsi_dot_fore_alpha(nidx,:) = [];
rslam(nidx) = []; %r for reduced
rng = 100;
cmap = turbo(rng);
dist = abs(rslam-E);
[sdist, sdidx] = sort(dist,'ascend');
ndist = (sdist - min(sdist))/(max(sdist) - min(sdist));
spPsi_dot_fore_alpha = pPsi_dot_fore_alpha(sdidx,:);
hold on
for i = 1:(length(rslam)-1)
    cidx = round(1 +(ndist(i)*(rng - 1)));
    plot(Sigma,spPsi_dot_fore_alpha(i,:),'Color',cmap(cidx,:),'LineWidth',1.5);
end
colorbar;
colormap(cmap);
clim([min(dist), max(dist)]);
ylim([-4, 4]);
set(ax1, 'FontSize', 16, 'TickLabelInterpreter', 'latex');
xlabel('$\sigma$','Interpreter','latex','FontSize',25)
ylabel('$\dot{c}_{E}(t)$','Interpreter','latex','FontSize',20)
ylabel(colorbar, '$|\bar{E}-E|$','Interpreter','latex','FontSize',15);
text(-0.1, 1, '(a)', 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 20,'Interpreter', 'latex')
box on


%%% === TILE 3
ax3 = nexttile(t,3);
hold(ax3, 'on');

rslam = slam;
nidx = find(slam>0); %if you want to remove E'<0 from plot
pPsi_dot_fore_beta = pPsi_dot_fore_full;
pPsi_dot_fore_beta(nidx,:) = [];
rslam(nidx) = []; %r for reduced
rng = 100;
cmap = turbo(rng);
dist = abs(rslam-E);
[sdist, sdidx] = sort(dist,'ascend');
ndist = (sdist - min(sdist))/(max(sdist) - min(sdist));
spPsi_dot_fore_beta = pPsi_dot_fore_beta(sdidx,:);
hold on
for i = 1:(length(rslam)-1)
    cidx = round(1 +(ndist(i)*(rng - 1)));
    plot(Sigma,spPsi_dot_fore_beta(i,:),'Color',cmap(cidx,:),'LineWidth',1.5);
end
colorbar;
colormap(cmap);
clim([min(sdist), max(sdist)]);
ylim([-4, 4]);
set(ax3, 'FontSize', 16, 'TickLabelInterpreter', 'latex');
xlabel('$\sigma$','Interpreter','latex','FontSize',25)
ylabel('$\dot{c}_{E}(t)$','Interpreter','latex','FontSize',20)
ylabel(colorbar, '$|\bar{E}-E|$','Interpreter','latex','FontSize',15);
text(-0.1, 1, '(c)', 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 20,'Interpreter', 'latex')
box on


%%% === TILE 2
ax2 = nexttile(t,2);
hold(ax2, 'on');
load('DS_Dyn1s_m1_E4_lam1_0_RG_k12_N1000.mat','E','sel2','slam','seve','Psi_dot_fore','Sigma')

pPsi_dot_fore = seve' * Psi_dot_fore;
pPsi_dot_fore_full = pPsi_dot_fore;
pPsi_dot_fore(sel2,:) = [];

rslam = slam;
nidx = find(slam<0); %if you want to remove E'<0 from plot
pPsi_dot_fore_alpha = pPsi_dot_fore_full;
pPsi_dot_fore_alpha(nidx,:) = [];
rslam(nidx) = []; %r for reduced
rng = 100;
cmap = turbo(rng);
dist = abs(rslam-E);
[sdist, sdidx] = sort(dist,'ascend');
ndist = (sdist - min(sdist))/(max(sdist) - min(sdist));
spPsi_dot_fore_alpha = pPsi_dot_fore_alpha(sdidx,:);
hold on
for i = 1:(length(rslam)-1)
    cidx = round(1 +(ndist(i)*(rng - 1)));
    plot(Sigma,spPsi_dot_fore_alpha(i,:),'Color',cmap(cidx,:),'LineWidth',1.5);
end
colorbar;
colormap(cmap);
clim([min(dist), max(dist)]);
ylim([-4, 4]);
set(ax2, 'FontSize', 16, 'TickLabelInterpreter', 'latex');
xlabel('$\sigma$','Interpreter','latex','FontSize',25)
ylabel('$\dot{c}_{E}(t)$','Interpreter','latex','FontSize',20)
ylabel(colorbar, '$|\bar{E}-E|$','Interpreter','latex','FontSize',15);
text(-0.1, 1, '(b)', 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 20,'Interpreter', 'latex')
box on


%%% === TILE 4
ax4 = nexttile(t,4);
hold(ax4, 'on');

rslam = slam;
nidx = find(slam>0); %if you want to remove E'<0 from plot
pPsi_dot_fore_beta = pPsi_dot_fore_full;
pPsi_dot_fore_beta(nidx,:) = [];
rslam(nidx) = []; %r for reduced
rng = 100;
cmap = turbo(rng);
dist = abs(rslam-E);
[sdist, sdidx] = sort(dist,'ascend');
ndist = (sdist - min(sdist))/(max(sdist) - min(sdist));
spPsi_dot_fore_beta = pPsi_dot_fore_beta(sdidx,:);
hold on
for i = 1:(length(rslam)-1)
    cidx = round(1 +(ndist(i)*(rng - 1)));
    plot(Sigma,spPsi_dot_fore_beta(i,:),'Color',cmap(cidx,:),'LineWidth',1.5);
end
colorbar;
colormap(cmap);
clim([min(dist), max(dist)]);
ylim([-4, 4]);
set(ax4, 'FontSize', 16, 'TickLabelInterpreter', 'latex');
xlabel('$\sigma$','Interpreter','latex','FontSize',25)
ylabel('$\dot{c}_{E}(t)$','Interpreter','latex','FontSize',20)
ylabel(colorbar, '$|\bar{E}-E|$','Interpreter','latex','FontSize',15);
text(-0.1, 1, '(d)', 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 20,'Interpreter', 'latex')
box on