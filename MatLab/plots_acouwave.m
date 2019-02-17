function plots_acouwave(varstr, ampl)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% plots_acouwave(varstr, ampl)
% Produces plots from RKLM low Mach fluid dynamics code hdf output
% for the paper "A semi-implicit numerical model for small-to-planetary scale atmospheric
% dynamics", acoustic wave test case with high (c=0.05*c_ref) or low (c=0.01*c_ref)
% initial amplitude ampl
%
%
% varstr    variable to plot
%           'p2_nodes', 'dp2_c', 'theta', 'dp2_nodes', 'dpdim', 'u',
%           'v', 'w', 'dY', 'rho', 'p', 'geopot', 'p2_c', 'rhoY'
%
% ampl      'high', 'low'
%
% Developed by R. Klein, FU Berlin, -2019
% Modified by T. Benacchio, Politecnico di Milano, 2019
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('/home/tommaso/work/code/matlab_packages/export_fig')
addpath('/home/tommaso/work/code/matlab_packages/Colormaps')

set(0,'DefaultFigureColor',[1 1 1])

linecolor      = 'default';
no_of_lines    = 25;
diff_rel_to_bottom = 0;

kmin = 0;
kmax = 2;
dk   = 1;

scalefactor = 1.0;
ncx = 256;
ncy = 10;
L   = 1.0 * scalefactor;  %
x0  = 0.5*L;
H   = 10.0;  %
aspect = [.1 1 1];
velosc = 1;  % velocity unit of RKLM code
showslice_hor = floor(ncy/2);

% auxiliary adjustments of grid parameters
dumsx = 2;
dumsy = 2;
transp    = 0;

folderstring = strcat('/home/tommaso/work/repos/RKLM_Reference/hdf_output/AcousticWave_', ampl);

folderstr = varstr;
ndummy = 2;
arraysize = [ncx ncy];

set(0,'defaulttextinterpreter','latex')
scrsz = get(0,'ScreenSize');

figure1 = figure('Position',[1 2*scrsz(4)/3 scrsz(4)/2 1*scrsz(4)/3]);

figure2 = figure('Position',[scrsz(4)/2 2*scrsz(4)/3 scrsz(4)/2 1*scrsz(4)/3]);
title(strcat('horizontal slice at j = ',num2str(showslice_hor)))

for k = kmin:dk:kmax
    kstr = num2str(k);
    filestr = strcat(folderstring,'/',folderstr,'/',varstr,'_00',kstr,'.hdf');
    v = hdfread(filestr, '/Data-Set-2', 'Index', {[1  1],[1  1],[arraysize(1)+dumsx*ndummy  arraysize(2)+dumsy*ndummy]});
    
    if transp == 1
        v = transpose(v);
    end
    
    velo = v;
    
    [nx, nz] = size(velo);
    nx = nx - 4;
    nz = nz - 4;
    
    dx = L/nx;
    dz = H/nz;
    
    x = linspace(x0 + 0.5*dx-(nx/2)*dx,x0 - 0.5*dx+(nx/2)*dx,nx);
    z = linspace(0.5*dz,-0.5*dz+nz*dz,nz);
    Yt = transpose(velo);
    th = Yt(3:1:nz+2, 3:1:nx+2);
    
    % Create filled contour
    figure(figure1)
    if (strcmp(varstr, 'u') || strcmp(varstr, 'v') || strcmp(varstr, 'w'))
        contourf(x,z,th*velosc,no_of_lines,'LineColor',linecolor);
    else
        if diff_rel_to_bottom
            th = th-th(3,:);
        end
        contourf(x,z,th,no_of_lines,'LineColor',linecolor);
    end
    colormap viridis
    colorbar('FontSize',14,'FontName','Helvetica')
    
    set(gca,'DataAspectRatio', aspect, 'FontSize',18,'FontName','Helvetica');
    axis tight;
    set(0,'defaulttextinterpreter','latex')
    
    
    figure(figure2)
    plot(th(showslice_hor,:), 'k', 'LineWidth', 2)
    
    xlabel('x [m]','FontSize',18,'Interpreter','latex');
    xlim([1 257])
    xticks([1 52 103 154 205 257])
    xticklabels({'0', '0.2', '0.4', '0.6', '0.8', '1'})
    if strcmp(varstr,'rho')
        ylabel('$\rho$ [$kg m^{-3}$]','FontSize',18,'Interpreter','latex');
        if strcmp(ampl,'high')
            ylim([0.8 1.2])
            yticks([0.9 1.0 1.1])
            yticklabels({'0.9', '1.0', '1.1'})
        elseif strcmp(ampl,'low')
            ylim([0.97 1.03])
            yticks([0.98 1.0 1.02])
            yticklabels({'0.98', '1.0', '1.02'})
        end
    elseif strcmp(varstr,'u')
        ylabel('$\rho u$ [$kg m^{-3}s^{-1}$]','FontSize',18,'Interpreter','latex');
        if strcmp(ampl,'high')
            ylim([-75 75])
            yticks([-50 0 50])
            yticklabels({'-50', '0', '50'})
        elseif strcmp(ampl,'low')
            ylim([-13 13])
            yticks([-10 0 10])
            yticklabels({'-10', '0', '10'})
        end
    end
    grid on
    fig=gcf;
    fig.Color = 'white';
    fig.InvertHardcopy = 'off';
    filename = sprintf('../RKLM_Reference/Doc/paper_2019/figures/AcousticWave_%s/%s/%s_snapshot%d.eps', ampl, varstr, varstr, k);
    print(filename, '-depsc')
    export_fig(filename, '-eps')
    
end