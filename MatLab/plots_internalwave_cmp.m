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
    contour_values_m = linspace(-0.00005, -0.00001,5);
    contour_values_p = linspace(0, 0.00015, 16);
else
    contour_values_m = linspace(-0.0002, -0.00004, 5)
    contour_values_p = linspace(0, 0.0002, 6);
end
[ccf1,hhf1]=contourf(x,z, th_C-th_PS, [min(min(th_C-th_PS)) contour_values_m, contour_values_p max(max(th_C-th_PS))], 'LineColor','k','LineWidth',1.0);
set(hhf1,'LineColor','none');
hold
contour(x,z, th_C-th_PS, [contour_values_p max(max(th_C-th_PS))],'LineColor','k','LineWidth',1.0)
[cc1,h1]=contour(x,z, th_C-th_PS, [min(min(th_C-th_PS)) contour_values_m], 'LineStyle', '--','LineColor','k','LineWidth',1.0);

% Take all the info from the contourline output argument:
i0 = 1;
i2 = 1;
while i0 <  length(cc1)
    i1 = i0+[1:cc1(2,i0)];
    zLevel(i2) = cc1(1,i0);
    hold on
    % And plot it with dashed lines:
    ph(i2) = plot(cc1(1,i1),cc1(2,i1),'k--','linewidth',1);
    i0 = i1(end)+1;
    i2 = i2+1;
end
% Scrap the contourlines:
delete(h1)

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
    contour_values_m=linspace(-0.00005, -0.00001, 5);
    contour_values_p=linspace(0, 0.00005, 6);
else
    contour_values_m=linspace(-0.000015, -0.000003, 5);
    contour_values_p=linspace(0.000003, 0.000015, 6);
end
[ccf1,hhf1]=contourf(x,z, th_C-th_HY, [min(min(th_C-th_HY)) contour_values_m contour_values_p max(max(th_C-th_HY))], 'LineColor','k','LineWidth',1.0);
set(hhf1,'LineColor','none');
hold
contour(x,z, th_C-th_HY, [contour_values_p max(max(th_C-th_HY))],'LineColor','k','LineWidth',1.0)
[cc1,h1]=contour(x,z, th_C-th_HY, [min(min(th_C-th_HY)) contour_values_m],'LineStyle', '--', 'LineColor','k','LineWidth',1.0);

% Take all the info from the contourline output argument:
i0 = 1;
i2 = 1;
while i0 <  length(cc1)
    i1 = i0+[1:cc1(2,i0)];
    zLevel(i2) = cc1(1,i0);
    hold on
    % And plot it with dashed lines:
    ph(i2) = plot(cc1(1,i1),cc1(2,i1),'k--','linewidth',1);
    i0 = i1(end)+1;
    i2 = i2+1;
end
% Scrap the contourlines:
delete(h1)


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