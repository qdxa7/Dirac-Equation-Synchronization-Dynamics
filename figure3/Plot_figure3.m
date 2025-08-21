tic

load('DS_Dyn1t_m1_E2_lam1_0_RG_k12_N1000_final');

pPsi = seve'*Psi_tot;
pPsi_dot = seve'*Psi_dot_tot;

%%Plot dynamical pPsi
[~,disp_dyn] = maxk(abs(pPsi_dot(:,length(time))),4);
disp_dyn = sort(disp_dyn);
dyn_pPsi = pPsi(disp_dyn,:);
slam_dyn = slam(disp_dyn);
dyn_pPsi_fig = figure;
dyn_pPsi_fig.Units = 'centimeters';
dyn_pPsi_fig.Position(3:4)=[35 20];
tiledlayout(2,2,'Padding','compact', 'TileSpacing','compact')
for i = 1:length(disp_dyn)
    nexttile
    if slam_dyn(i)==E
        ylab = "$\sin(c_{\overline{E}})$";
        colour = 'r';
    else
        ylab = "$\sin(c_{E})$";
        colour = 'b';
    end
    label = {'a' 'b' 'c' 'd'};
    plot(time, sin(dyn_pPsi(i,:)),colour,'LineWidth',3); %,'DisplayName',sprintf(leg,slam_dyn(i)))
    set(gca,'TickLabelInterpreter', 'latex','FontWeight','normal','FontSize',16);
    xlabel('$t$','Interpreter', 'latex','FontSize',25)
    ylabel(ylab,'Interpreter', 'latex','FontSize',25)
    ylim([-1 1])
    text(-0.12, 1, sprintf('(%c)',label{i}), 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 20,'Interpreter', 'latex')
end

toc