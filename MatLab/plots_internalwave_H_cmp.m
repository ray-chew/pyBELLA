%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% plots_internalwave_H_cmp
% Produces plots with difference between configurations for the hydrostatic version of the
% hydrostatic SK94 internal gravity wave from RKLM low Mach fluid dynamics code hdf output
% for the paper "A semi-implicit numerical model for small-to-planetary scale atmospheric
% dynamics", wave test case
%
%
% Developed by R. Klein, FU Berlin, -2019
% Modified by T. Benacchio, Politecnico di Milano, 2019
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

addpath('./export_fig')

[x,z,th_C] = plots_internalwave('dY', 'H');
[x,z,th_PS] = plots_internalwave('dY', 'H_psinc');
[x,z,th_HY] = plots_internalwave('dY', 'H_hyd');
figure(1)

set(0,'DefaultFigureColor',[1 1 1])

L   = 6000.0;  %
H   = 10.0;  %

aspect = [L/H/3 1 1];

set(0,'defaulttextinterpreter','latex')
scrsz = get(0,'ScreenSize');

figure1 = figure('Position',[1 2*scrsz(4)/3 scrsz(4)/2 1*scrsz(4)/3]);

figure(figure1)
contour_values=linspace(-0.00005, 0.00015, 21);
contourf(x,z, th_C-th_PS, [min(min(th_C-th_PS)) contour_values max(max(th_C-th_PS))], 'LineColor','k','LineWidth',1.0)
hold
contour(x,z, th_C-th_PS, [min(min(th_C-th_PS)) contour_values max(max(th_C-th_PS))],'LineColor','k','LineWidth',1.0)

set(gca,'DataAspectRatio', aspect, 'FontSize',18,'FontName','Helvetica');
axis tight;

set(0,'defaulttextinterpreter','latex')
set(0,'DefaultFigureColor',[1 1 1])
xlabel('x [km]','FontSize',18,'Interpreter','latex');
ylabel('z [km]','FontSize',18,'Interpreter','latex');
yticks([0 2 4 6 8 10])
colormap viridis

xlim([-3000 3000])
xticks([-3000 -2000 -1000 0 1000 2000 3000])
xticklabels({'0', '1000', '2000', '3000', '4000', '5000', '6000'})
colorbar('FontSize',14,'FontName','Helvetica');

fig=gcf;
fig.Color = 'white';
fig.InvertHardcopy = 'off';
filename = sprintf('../RKLM_Reference/Doc/paper_2019/figures/InternalWave_H/dY/CMP_vs_PS.eps');
print(filename, '-depsc')
export_fig(filename, '-eps')

figure2 = figure('Position',[1 2*scrsz(4)/3 scrsz(4)/2 1*scrsz(4)/3]);

figure(figure2)
contour_values=linspace(-0.00005, 0.00005, 11);
contourf(x,z, th_C-th_HY, [min(min(th_C-th_HY)) contour_values max(max(th_C-th_HY))], 'LineColor','k','LineWidth',1.0)
hold
contour(x,z, th_C-th_HY, [min(min(th_C-th_HY)) contour_values max(max(th_C-th_HY))],'LineColor','k','LineWidth',1.0)

set(gca,'DataAspectRatio', aspect, 'FontSize',18,'FontName','Helvetica');
axis tight;

set(0,'defaulttextinterpreter','latex')
set(0,'DefaultFigureColor',[1 1 1])
xlabel('x [km]','FontSize',18,'Interpreter','latex');
ylabel('z [km]','FontSize',18,'Interpreter','latex');
yticks([0 2 4 6 8 10])
colormap viridis

xlim([-3000 3000])
xticks([-3000 -2000 -1000 0 1000 2000 3000])
xticklabels({'0', '1000', '2000', '3000', '4000', '5000', '6000'})
colorbar('FontSize',14,'FontName','Helvetica');

fig=gcf;
fig.Color = 'white';
fig.InvertHardcopy = 'off';
filename = sprintf('../RKLM_Reference/Doc/paper_2019/figures/InternalWave_H/dY/CMP_vs_HY.eps');
print(filename, '-depsc')
export_fig(filename, '-eps')
