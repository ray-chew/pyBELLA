function plots_vortex(varstr, resol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% plots_vortex(varstr, ext)
% Produces plots from RKLM low Mach fluid dynamics code hdf output
% for the paper "A semi-implicit numerical model for small-to-planetary scale atmospheric
% dynamics", rotating Gresho vortex case.
%
%
% varstr    variable to plot
%           'p2_nodes', 'dp2_c', 'theta', 'dp2_nodes', 'dpdim', 'u',
%           'v', 'w', 'dY', 'rho', 'p', 'geopot', 'p2_c', 'rhoY'
%
% resol     resolution in number of points, eg 96, 192
%
% Developed by R. Klein, FU Berlin, -2019
% Modified by T. Benacchio, Politecnico di Milano, 2019
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('/home/tommaso/work/code/matlab_packages/export_fig')
addpath('/home/tommaso/work/code/matlab_packages/Colormaps')

set(0,'DefaultFigureColor',[1 1 1])

linecolor      = 'default';  % 'k', 'default' ...

kmin = 0;
kmax = 3;
dk   = 1;


ncx = resol;
ncy = resol;
L   = 4.0;
x0  = 0.5;
H   = 1.0;
aspect = [2 2 2];
velosc = 1;
showslice_hor = floor(ncy/2);
showslice_ver = floor(ncx/2);

if strcmp(varstr, 'rho')
    contour_values = linspace(0.525, 0.975, 19);
elseif strcmp(varstr, 'p2_n')
    contour_values = linspace(-0.001, -.005, 10);
    nnx=ncx+1;
    nny=ncy+1;
end    
% auxiliary adjustments of grid parameters
dumsx = 1;
dumsy = 2;

folderstring = strcat('/home/tommaso/work/repos/RKLM_Reference/hdf_output/TravellingVortex_3D_', num2str(resol));

% cell-centered fields
folderstr = varstr;
if strcmp(varstr, 'p2_n')
    folderstr='p2_nodes';
end
ndummy = 2;
if strcmp(varstr, 'rho')
    arraysize = [ncx ncy];
elseif strcmp(varstr, 'p2_n')
    arraysize = [nnx nny];
end

set(0,'defaulttextinterpreter','latex')
scrsz = get(0,'ScreenSize');

figure1 = figure('Position',[1 2*scrsz(4)/3 scrsz(4)/2 1*scrsz(4)/3]);

for k = kmin:dk:kmax
    kstr = num2str(k);
    filestr = strcat(folderstring,'/',folderstr,'/',varstr,'_00',kstr,'.hdf');
    
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
    
    % Create filled contour
    figure(figure1)
    contourf(x,z,th,[min(min(th)) contour_values max(max(th))],'LineColor',linecolor);
    hold on
    contour(x,z,th,[min(min(th)) contour_values max(max(th))],'LineColor','k');
    
    colormap viridis
    colorbar('FontSize',14,'FontName','Helvetica')

    set(gca,'DataAspectRatio', aspect, 'FontSize',18,'FontName','Helvetica');
    axis tight;
    
    set(0,'defaulttextinterpreter','latex')
    set(0,'DefaultFigureColor',[1 1 1])
    xlim([-1.5 2.5])
    xticks([-1.5 0.5 2.5])
    xticklabels([{'0', '0.5', '1'}'])
    xlabel('x [m]','FontSize',18,'Interpreter','latex');
    ylabel('z [m]','FontSize',18,'Interpreter','latex');
    yticks([0 0.5 1])
    axis square
    
    fig=gcf;
    fig.Color = 'white';
    fig.InvertHardcopy = 'off';
    filename = sprintf('../RKLM_Reference/Doc/paper_2019/figures/TravellingVortex_%s/%s/%s_snapshot%d.eps', num2str(resol), varstr, varstr, k);
    print(filename, '-depsc')
    export_fig(filename, '-eps')
    
end