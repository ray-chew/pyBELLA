function [x,z,th]=plots_internalwave(varstr, ext)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% plots_internalwave(varstr, ext)
% Produces plots from RKLM low Mach fluid dynamics code hdf output
% for the paper "A semi-implicit numerical model for small-to-planetary scale atmospheric
% dynamics", wave test case with nonhydrostatic (NH), hydrostatic(H), or
% planetary(P) extension ext
%
%
% varstr    variable to plot
%           'p2_nodes', 'dp2_c', 'theta', 'dp2_nodes', 'dpdim', 'u',
%           'v', 'w', 'dY', 'rho', 'p', 'geopot', 'p2_c', 'rhoY'
%
% ext      Extension of the domain
%           'NH', 'H', 'H_hyd', 'H_psinc', 'P', 'P_hyd', 'P_psinc'
%
% Developed by R. Klein, FU Berlin, -2019
% Modified by T. Benacchio, Politecnico di Milano, 2019
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('./export_fig')

set(0,'DefaultFigureColor',[1 1 1])

dtheta = 0.5e-3;
contour_values = [-5*dtheta, -4*dtheta, -3*dtheta, -2*dtheta, -dtheta, 0.0, dtheta, 2*dtheta, 3*dtheta, 4*dtheta, 5*dtheta];

kmin = 0;
kmax = 2;
dk   = 1;

if strcmp(ext, 'NH')
    scalefactor = 1.0;  % Skamarock-Klemp-1994 Fig.1
elseif strcmp(ext, 'H') || strcmp(ext, 'H_psinc') || strcmp(ext, 'H_hyd')
    scalefactor = 20.0;  % Skamarock-Klemp-1994 Fig.3
elseif strcmp(ext, 'P') || strcmp(ext, 'P_psinc') || strcmp(ext, 'P_hyd')
    scalefactor = 160.0;   % new, very long wave test
end


ncx = 301;
ncy = 10;
L   = 300.0 * scalefactor;  %
x0  = 0.0;
H   = 10.0;  %
aspect = [L/H/3 1 1];

% auxiliary adjustments of grid parameters
dumsx = 2;
dumsy = 2;

folderstring = strcat('../hdf_output/InternalWave_', ext);

folderstr = varstr;
ndummy = 2;
arraysize = [ncx ncy];


set(0,'defaulttextinterpreter','latex')
scrsz = get(0,'ScreenSize');

figure1 = figure('Position',[1 2*scrsz(4)/3 scrsz(4)/2 1*scrsz(4)/3]);

for k = kmin:dk:kmax
    kstr = num2str(k);
    filestr = strcat(folderstring,'/',folderstr,'/',varstr,'_00',kstr,'.hdf');
    v = hdfread(filestr, '/Data-Set-2', 'Index', {[1  1],[1  1],[arraysize(1)+dumsx*ndummy  arraysize(2)+dumsy*ndummy]});
    
    velo = v;
    
    [nx, nz] = size(velo);
    nx = nx - 4;
    nz = nz - 4;
    
    dx = L/nx;
    dz = H/nz;
    
    x = linspace(x0 + 0.5*dx-(nx/2)*dx,x0 - 0.5*dx+(nx/2)*dx,nx);
    z = linspace(0.5*dz,-0.5*dz+nz*dz,nz);
    Yt = transpose(velo);
    th = 300*Yt(3:1:nz+2, 3:1:nx+2);
    
    % Create unfilled contour
    figure(figure1)
    if k==0
        contour_values = linspace(0, 0.01, 11);
    else
        if strcmp(ext, 'NH') || strcmp(ext, 'H') || strcmp(ext, 'H_psinc') || strcmp(ext, 'H_hyd')
           contour_values = [-5*dtheta, -4*dtheta, -3*dtheta, -2*dtheta, -dtheta, 0.0, dtheta, 2*dtheta, 3*dtheta, 4*dtheta, 5*dtheta];
        else
           contour_values = 2.0*[-5*dtheta, -4*dtheta, -3*dtheta, -2*dtheta, -dtheta, 0.0, dtheta, 2*dtheta, 3*dtheta, 4*dtheta, 5*dtheta];
        end       
    end
    contourf(x,z,th,[min(min(th)) contour_values max(max(th))],'LineColor','k','LineWidth',1.0);
    %contour(x,z,max(0.0,th),contour_values,'LineColor','k','LineWidth',1.0);
    hold
    contour(x,z,th,[min(min(th)) contour_values max(max(th))],'LineColor','k','LineWidth',1.0);
%    contour(x,z,min(0.0,th),contour_values,'LineColor','k');
    hold
    
    set(gca,'DataAspectRatio', aspect, 'FontSize',14,'FontName','Helvetica');
    axis tight;
    
    set(0,'defaulttextinterpreter','latex')
    set(0,'DefaultFigureColor',[1 1 1])
    xlabel('x [km]','FontSize',18,'Interpreter','latex');
    ylabel('z [km]','FontSize',18,'Interpreter','latex');
    yticks([0 2 4 6 8 10])
    colormap viridis
    
    if k==0
        colorbar('FontSize',14,'FontName','Helvetica')
    end
    
    switch ext
        case 'NH'
            xlim([-150 150])
            xticks([-150 -100 -50 0 50 100 150])
            xticklabels({'0', '50', '100', '150', '200', '250', '300'})
            if k~=0; colorbar('FontSize',14,'FontName','Helvetica'); end
        case 'H'
            xlim([-3000 3000])
            xticks([-3000 -2000 -1000 0 1000 2000 3000])
            xticklabels({'0', '1000', '2000', '3000', '4000', '5000', '6000'})
            if k~=0; colorbar('FontSize',14,'FontName','Helvetica'); end
        case 'H_psinc'
            xlim([-3000 3000])
            xticks([-3000 -2000 -1000 0 1000 2000 3000])
            xticklabels({'0', '1000', '2000', '3000', '4000', '5000', '6000'})
            if k~=0; colorbar('FontSize',14,'FontName','Helvetica'); end
        case 'H_hyd'
            xlim([-3000 3000])
            xticks([-3000 -2000 -1000 0 1000 2000 3000])
            xticklabels({'0', '1000', '2000', '3000', '4000', '5000', '6000'})
            if k~=0; colorbar('FontSize',14,'FontName','Helvetica'); end
        case 'P'
            xlabel('x [$10^3$ km]','FontSize',18,'Interpreter','latex');
            xlim([-24000 24000])
            xticks([-24000 -16000 -8000 0 8000 16000 24000])
            xticklabels({'0', '8', '16', '24', '32', '40', '48'})
            if k~=0; colorbar('FontSize',14,'FontName','Helvetica'); end
        case 'P_psinc'
            xlabel('x [$10^3$ km]','FontSize',18,'Interpreter','latex');
            xlim([-24000 24000])
            xticks([-24000 -16000 -8000 0 8000 16000 24000])
            xticklabels({'0', '8', '16', '24', '32', '40', '48'})
            if k~=0; colorbar('FontSize',14,'FontName','Helvetica'); end
        case 'P_hyd'
            xlabel('x [$10^3$ km]','FontSize',18,'Interpreter','latex');
            xlim([-24000 24000])
            xticks([-24000 -16000 -8000 0 8000 16000 24000])
            xticklabels({'0', '8', '16', '24', '32', '40', '48'})
            if k~=0; colorbar('FontSize',14,'FontName','Helvetica'); end
        case default
            disp('Error: incorrect test_case input, see help make_plots. Exiting.')
            exit;
    end
    
    fig=gcf;
    fig.Color = 'white';
    fig.InvertHardcopy = 'off';
    filename = sprintf('../RKLM_Reference/Doc/paper_2019/figures/InternalWave_%s/%s/%s_snapshot%d.eps', ext, varstr, varstr, k);
    print(filename, '-depsc')
    export_fig(filename, '-eps')
    
end