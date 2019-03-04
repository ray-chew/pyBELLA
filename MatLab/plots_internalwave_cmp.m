function plots_internalwave_cmp(conf)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% plots_internalwave_cmp
% Produces plots with difference between configurations for the hydrostatic version of the
% hydrostatic SK94 internal gravity wave from RKLM low Mach fluid dynamics code hdf output
% for the paper "A semi-implicit numerical model for small-to-planetary scale atmospheric
% dynamics", wave test case
%
% conf configuration at which to plot the differences
%  
% 'H', 'P'
%
% Developed by R. Klein, FU Berlin, -2019
% Modified by T. Benacchio, Politecnico di Milano, 2019
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('./export_fig')

[x,z,th_C] = plots_internalwave('dY', conf);
[x,z,th_PS] = plots_internalwave('dY', strcat(conf, '_psinc'));
[x,z,th_HY] = plots_internalwave('dY', strcat(conf, '_hyd'));
figure(1)

set(0,'DefaultFigureColor',[1 1 1])

if strcmp(conf, 'H')
    L=6000.0;
else
    L=48000.0;
end
H=10.0;

aspect = [L/H/3 1 1];

set(0,'defaulttextinterpreter','latex')
scrsz = get(0,'ScreenSize');

figure1 = figure('Position',[1 2*scrsz(4)/3 scrsz(4)/2 1*scrsz(4)/3]);

figure(figure1)
if strcmp(conf, 'H')
    contour_values=linspace(-0.00005, 0.00015, 21);
else
    contour_values=linspace(-0.0002, 0.0002, 11);
end
contourf(x,z, th_C-th_PS, [min(min(th_C-th_PS)) contour_values max(max(th_C-th_PS))], 'LineColor','k','LineWidth',1.0)
%contourf(x,z, th_C-th_PS, 21, 'LineColor','k','LineWidth',1.0)
hold
contour(x,z, th_C-th_PS, [min(min(th_C-th_PS)) contour_values max(max(th_C-th_PS))],'LineColor','k','LineWidth',1.0)
%[C,h]=contour(x,z, th_C-th_PS, 21,'LineColor','k','LineWidth',1.0)
%clabel(C,h)

set(gca,'DataAspectRatio', aspect, 'FontSize',14,'FontName','Helvetica');
axis tight;

set(0,'defaulttextinterpreter','latex')
set(0,'DefaultFigureColor',[1 1 1])
ylabel('z [km]','FontSize',18,'Interpreter','latex');
yticks([0 2 4 6 8 10])
colormap viridis

if strcmp(conf, 'H')
    xlim([-3000 3000])
    xlabel('x [km]','FontSize',18,'Interpreter','latex');
    xticks([-3000 -2000 -1000 0 1000 2000 3000])
    xticklabels({'0', '1000', '2000', '3000', '4000', '5000', '6000'})
else
    xlim([-24000 24000])
    xlabel('x [$10^3$ km]','FontSize',18,'Interpreter','latex');
    xticks([-24000 -16000 -8000 0 8000 16000 24000])
    xticklabels({'0', '8', '16', '24', '32', '40', '48'})
end
colorbar('FontSize',14,'FontName','Helvetica');

fig=gcf;
fig.Color = 'white';
fig.InvertHardcopy = 'off';
filename = sprintf('../RKLM_Reference/Doc/paper_2019/figures/InternalWave_%s/dY/CMP_vs_PS.eps', conf);
print(filename, '-depsc')
export_fig(filename, '-eps')

figure2 = figure('Position',[1 2*scrsz(4)/3 scrsz(4)/2 1*scrsz(4)/3]);

figure(figure2)
if strcmp(conf, 'H')
    contour_values=linspace(-0.00005, 0.00005, 11);
else
    contour_values=linspace(-0.000015, 0.000015, 11);
end
contourf(x,z, th_C-th_HY, [min(min(th_C-th_HY)) contour_values max(max(th_C-th_HY))], 'LineColor','k','LineWidth',1.0)
%contourf(x,z, th_C-th_HY, 11, 'LineColor','k','LineWidth',1.0)
hold
contour(x,z, th_C-th_HY, [min(min(th_C-th_HY)) contour_values max(max(th_C-th_HY))],'LineColor','k','LineWidth',1.0)
%contour(x,z, th_C-th_HY, 11,'LineColor','k','LineWidth',1.0)

set(gca,'DataAspectRatio', aspect, 'FontSize',14,'FontName','Helvetica');
axis tight;

set(0,'defaulttextinterpreter','latex')
set(0,'DefaultFigureColor',[1 1 1])
xlabel('x [km]','FontSize',18,'Interpreter','latex');
ylabel('z [km]','FontSize',18,'Interpreter','latex');
yticks([0 2 4 6 8 10])
colormap viridis

if strcmp(conf, 'H')
    xlim([-3000 3000])
    xlabel('x [km]','FontSize',18,'Interpreter','latex');
    xticks([-3000 -2000 -1000 0 1000 2000 3000])
    xticklabels({'0', '1000', '2000', '3000', '4000', '5000', '6000'})
else
    xlim([-24000 24000])
    xlabel('x [$10^3$ km]','FontSize',18,'Interpreter','latex');
    xticks([-24000 -16000 -8000 0 8000 16000 24000])
    xticklabels({'0', '8', '16', '24', '32', '40', '48'})
end
colorbar('FontSize',14,'FontName','Helvetica');

fig=gcf;
fig.Color = 'white';
fig.InvertHardcopy = 'off';
filename = sprintf('../RKLM_Reference/Doc/paper_2019/figures/InternalWave_%s/dY/CMP_vs_HY.eps', conf);
print(filename, '-depsc')
export_fig(filename, '-eps')
end