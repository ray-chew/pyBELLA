function plots_vortex_convergence(varstr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% plots_vortex_convergence(varstr)
% Produces convergence plots from RKLM low Mach fluid dynamics code hdf output
% for the paper "A semi-implicit numerical model for small-to-planetary scale atmospheric
% dynamics", rotating Gresho vortex case.
%
%
% varstr    variable to plot
%           'p2_nodes', 'dp2_c', 'theta', 'dp2_nodes', 'dpdim', 'u',
%           'v', 'w', 'dY', 'rho', 'p', 'geopot', 'p2_c', 'rhoY'
%
%
% Developed by R. Klein, FU Berlin, -2019
% Modified by T. Benacchio, Politecnico di Milano, 2019
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('/export_fig')

set(0,'DefaultFigureColor',[1 1 1])

linecolor      = 'default';  % 'k', 'default' ...

kmin = 0;
kmax = 1;

x = [48 96 192 384];



ncx = 192;
ncy = 192;

elseif strcmp(varstr, 'p2_n')
    nnx=ncx+1;
    nny=ncy+1;
end

% auxiliary adjustments of grid parameters
dumsx = 1;
dumsy = 2;

L   = 4.0;
x0  = 0.5;
H   = 1.0;
ndummy = 2;

% 48x48

ncx = 48;
ncy = 48;
nnx = ncy+1;
nny = ncy+1;

if strcmp(varstr, 'rho')
    arraysize = [ncx ncy];
elseif strcmp(varstr, 'p2_n')
    arraysize = [nnx nny];
end

folderstring = strcat('/home/tommaso/work/repos/RKLM_Reference/hdf_output/TravellingVortex_3D_48');
filestr = strcat(folderstring,'/',folderstr,'/',varstr,'_001.hdf');

v = hdfread(filestr, '/Data-Set-2', 'Index', {[1  1],[1  1],[arraysize(1)+dumsx*ndummy  arraysize(2)+dumsy*ndummy]});

velo=v;
[nx, nz] = size(velo);
nx = nx - 4;
nz = nz - 4;

Yt = transpose(velo);
% Copy-pasting the first and last rows to get a square matrix
th1(1, :) = Yt(3, 3:1:nx+2);
th1(2:nz+1, :) = Yt(3:1:nz+2, 3:1:nx+2);
th1(nz+2, :) = Yt(nz+2, 3:1:nx+2);

% 96x96
ncx = 2*ncx;
ncy = 2*ncy;
nnx = ncy+1;
nny = ncy+1;

if strcmp(varstr, 'rho')
    arraysize = [ncx ncy];
elseif strcmp(varstr, 'p2_n')
    arraysize = [nnx nny];
end


folderstring = strcat('/home/tommaso/work/repos/RKLM_Reference/hdf_output/TravellingVortex_3D_96');
filestr = strcat(folderstring,'/',folderstr,'/',varstr,'_001.hdf');

v = hdfread(filestr, '/Data-Set-2', 'Index', {[1  1],[1  1],[arraysize(1)+dumsx*ndummy  arraysize(2)+dumsy*ndummy]});

velo=v;
[nx, nz] = size(velo);
nx = nx - 4;
nz = nz - 4;

Yt = transpose(velo);
% Copy-pasting the first and last rows to get a square matrix
th2(1, :) = Yt(3, 3:1:nx+2);
th2(2:nz+1, :) = Yt(3:1:nz+2, 3:1:nx+2);
th2(nz+2, :) = Yt(nz+2, 3:1:nx+2);

size(th2)

% 192x192
ncx = 2*ncx;
ncy = 2*ncy;
nnx = ncy+1;
nny = ncy+1;

if strcmp(varstr, 'rho')
    arraysize = [ncx ncy];
elseif strcmp(varstr, 'p2_n')
    arraysize = [nnx nny];
end
folderstring = strcat('/home/tommaso/work/repos/RKLM_Reference/hdf_output/TravellingVortex_3D_192');
filestr = strcat(folderstring,'/',folderstr,'/',varstr,'_001.hdf');

v = hdfread(filestr, '/Data-Set-2', 'Index', {[1  1],[1  1],[arraysize(1)+dumsx*ndummy  arraysize(2)+dumsy*ndummy]});

velo=v;
[nx, nz] = size(velo);
nx = nx - 4;
nz = nz - 4;

% Copy-pasting the first and last rows to get a square matrix
Yt = transpose(velo);
th3(1, :) = Yt(3, 3:1:nx+2);
th3(2:nz+1, :) = Yt(3:1:nz+2, 3:1:nx+2);
th3(nz+2, :) = Yt(nz+2, 3:1:nx+2);

size(th3)


% 384x384
ncx = 2*ncx;
ncy = 2*ncy;
nnx = ncy+1;
nny = ncy+1;

if strcmp(varstr, 'rho')
    arraysize = [ncx ncy];
elseif strcmp(varstr, 'p2_n')
    arraysize = [nnx nny];
end

folderstring = strcat('/home/tommaso/work/repos/RKLM_Reference/hdf_output/TravellingVortex_3D_384');
filestr = strcat(folderstring,'/',folderstr,'/',varstr,'_001.hdf');

v = hdfread(filestr, '/Data-Set-2', 'Index', {[1  1],[1  1],[arraysize(1)+dumsx*ndummy  arraysize(2)+dumsy*ndummy]});

velo=v;
[nx, nz] = size(velo);
nx = nx - 4;
nz = nz - 4;

Yt = transpose(velo);
% Copy-pasting the first and last rows to get a square matrix
th4(1, :) = Yt(3, 3:1:nx+2);
th4(2:nz+1, :) = Yt(3:1:nz+2, 3:1:nx+2);
th4(nz+2, :) = Yt(nz+2, 3:1:nx+2);

size(th4)

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

    filestr = strcat(folderstring,'/',folderstr,'/',varstr,'_001.hdf');
    
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
    
    size(th);
    
    % Create filled contour
    figure(figure1)
    
    fig=gcf;
    fig.Color = 'white';
    fig.InvertHardcopy = 'off';
    filename = sprintf('../RKLM_Reference/Doc/paper_2019/figures/TravellingVortex/%s/%s_snapshot%d.eps', varstr, varstr, k);
    print(filename, '-depsc')
    export_fig(filename, '-eps')
    