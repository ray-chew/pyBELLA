function plots_straka(varstr, resol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% plots_straka(varstr, ext)
% Produces plots from RKLM low Mach fluid dynamics code hdf output
% for the paper "A semi-implicit numerical model for small-to-planetary scale atmospheric
% dynamics",density current test case
%
%
% varstr    variable to plot
%           'p2_nodes', 'dp2_c', 'theta', 'dp2_nodes', 'dpdim', 'u',
%           'v', 'w', 'dY', 'rho', 'p', 'geopot', 'p2_c', 'rhoY'
%
% resol     number of points in the horizontal direction, eg 1025 for 50m
%           resolution
%
% Developed by R. Klein, FU Berlin, -2019
% Modified by T. Benacchio, Politecnico di Milano, 2019
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('./export_fig')

set(0,'DefaultFigureColor',[1 1 1])

linecolor      = 'default';  % 'k', 'default' ...

kmin = 0;
kmax = 3;
dk   = 1;


ncx = 51200/resol + 1;
ncy = 6400/resol;
L  = 51.2;  
x0 = 0.0*L;
H  = 6.4;  
aspect = [1 1 1];
dtheta = 1.0;
contour_values = linspace(-16.5*dtheta,-0.5*dtheta,17);

% auxiliary adjustments of grid parameters
dumsx = 2;
dumsy = 2;

folderstring = strcat('../hdf_output/Straka_3D_', num2str(resol), 'm');

% cell-centered fields
folderstr = varstr;
ndummy = 2;
arraysize = [ncx ncy];

set(0,'defaulttextinterpreter','latex')
scrsz = get(0,'ScreenSize');

figure1 = figure('Position',[1 2*scrsz(4)/3 scrsz(4)/2 1*scrsz(4)/3]);

for k = kmin:dk:kmax
    kstr = num2str(k);
    filestr = strcat(folderstring,'/',folderstr,'/',varstr,'_00',kstr,'.hdf');
    filestr
    v = hdfread(filestr, '/Data-Set-2', 'Index', {[1  1],[1  1],[arraysize(1)+dumsx*ndummy  arraysize(2)+dumsy*ndummy]});
    
    velo=v;
    [nx, nz] = size(velo);
    nx = nx - 4;
    nz = nz - 4;
    
    dx = L/nx;
    dz = H/nz;
    
    x = linspace(x0 + 0.5*dx-(nx/2)*dx,x0 - 0.5*dx+(nx/2)*dx,nx);
    z = linspace(0.5*dz,-0.5*dz+nz*dz,nz);
    Yt = transpose(velo);
    th = Yt(3:1:nz+2, 3:1:nx+2);
    
    if k==kmax, fprintf('min(theta'')=%2.4f K\nmax(theta'')=%2.4f K\n', min(min(th))*300, max(max(th))*300), end
    
    % Create filled contour
    figure(figure1)
    contourf(x,z,300*th,contour_values,'LineColor',linecolor);
    hold on
    contour(x,z,300*th,contour_values,'LineColor','k','LineWidth',1.0);
    colormap viridis
    colorbar('FontSize',14,'FontName','Helvetica');

    set(gca,'DataAspectRatio', aspect, 'FontSize',18,'FontName','Helvetica');
    axis tight;
    
    set(0,'defaulttextinterpreter','latex')
    set(0,'DefaultFigureColor',[1 1 1])
    xlim([0 19.2])
    xticks([0 2 4 6 8 10 12 14 16 18])
    xlabel('x [km]','FontSize',18,'Interpreter','latex');
    ylabel('z [km]','FontSize',18,'Interpreter','latex');
    ylim([min(z) 5])
    yticks([0 1 2 3 4])
    
    fig=gcf;
    fig.Color = 'white';
    fig.InvertHardcopy = 'off';
    filename = sprintf('../RKLM_Reference/Doc/paper_2019/figures/Straka_%sm/%s/%s_snapshot%d.eps', num2str(resol), varstr, varstr, k);
    print(filename, '-depsc')
    export_fig(filename, '-eps')

    if k==kmax
       hold off 
       figure(k)
       [cc,handle]=contour(x,z,300*th,[-1 -1],'LineColor','k');
       fprintf('Front position:%5d m\n', floor(1000*max(cc(1,:))));
    end

    
end