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

% auxiliary adjustments of grid parameters
dumsx = 1;
dumsy = 2;

L   = 4.0;
x0  = 0.5;
H   = 1.0;
ndummy = 2;

% cell-centered fields
folderstr = varstr;
if strcmp(varstr, 'p2_n')
    folderstr='p2_nodes';
end


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
th1(:, 1)      = Yt(3:1:nz+2, 3);
th1(:, 2:nx+1) = Yt(3:1:nz+2, 3:1:nx+2);
th1(:, nx+2)   = Yt(3:1:nz+2, nx+2);

% Chop off first row and column for simplicity for nodal variables
if strcmp(varstr, 'p2_n')
  th1=th1(2:nz, 2:nx+2);    
end

size(th1)

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
th2(:, 1)      = Yt(3:1:nz+2, 3);
th2(:, 2:nx+1) = Yt(3:1:nz+2, 3:1:nx+2);
th2(:, nx+2)  = Yt(3:1:nz+2, nx+2);

% Chop off first row and column for simplicity for nodal variables
if strcmp(varstr, 'p2_n')
  th2=th2(2:nz, 2:nx+2);    
end


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
th3(:, 1)      = Yt(3:1:nz+2, 3);
th3(:, 2:nx+1) = Yt(3:1:nz+2, 3:1:nx+2);
th3(:, nx+2)  = Yt(3:1:nz+2, nx+2);

% Chop off first row and column for simplicity for nodal variables
if strcmp(varstr, 'p2_n')
  th3=th3(2:nz, 2:nx+2);    
end


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
th4(:, 1)      = Yt(3:1:nz+2, 3);
th4(:, 2:nx+1) = Yt(3:1:nz+2, 3:1:nx+2);
th4(:, nx+2)  = Yt(3:1:nz+2, nx+2);

% Chop off first row and column for simplicity for nodal variables
if strcmp(varstr, 'p2_n')
  th4=th4(2:nz, 2:nx+2);    
end

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
th5(:, 1)      = Yt(3:1:nz+2, 3);
th5(:, 2:nx+1) = Yt(3:1:nz+2, 3:1:nx+2);
th5(:, nx+2)  = Yt(3:1:nz+2, nx+2);

% Chop off first row and column for simplicity for nodal variables
if strcmp(varstr, 'p2_n')
  th5=th5(2:nz, 2:nx+2);    
end

% ========================================================
% Creating coarse-grid versions of fine-grid solution
% ========================================================


for j = 1:2:size(th4, 2)-1
    for i = 1:2:size(th5, 1)-1
        th5_av4(floor(i/2)+1, floor(j/2)+1) = .25*(th5(i, j)+ th5(i, j+1)+ th5(i+1, j)+th5(i+1, j+1));
    end
end

for j = 1:4:size(th5, 2)-1
    for i = 1:4:size(th5, 1)-1
        th5_av3(floor(i/4)+1, floor(j/4)+1) = .25*(th5(i, j) + th5(i, j+3)+ th5(i+3, j)+th5(i+3, j+3));
    end
end

for j = 1:8:size(th5, 2)-1
    for i = 1:8:size(th5, 1)-1
        th5_av2(floor(i/8)+1, floor(j/8)+1) = .25*(th5(i, j) + th5(i, j+7) + th5(i+7, j) + th5(i+7, j+7));
    end
end

for j = 1:16:size(th5, 2)-1
    for i = 1:16:size(th5, 1)-1
        th5_av1(floor(i/16)+1, floor(j/16)+1) = .25*(th5(i, j) + th5(i, j+15) + th5(i+15, j) + th5(i+15, j+15));
    end
end

rel_L2_err_1_th = norm(th1-th5_av1, 2)/norm(th1, 2);
rel_L2_err_2_th = norm(th2-th5_av2, 2)/norm(th2, 2);
rel_L2_err_3_th = norm(th3-th5_av3, 2)/norm(th3, 2);
rel_L2_err_4_th = norm(th4-th5_av4, 2)/norm(th4, 2);

rel_Linf_err_1_th = norm(th1-th5_av1, inf)/norm(th1, inf);
rel_Linf_err_2_th = norm(th2-th5_av2, inf)/norm(th2, inf);
rel_Linf_err_3_th = norm(th3-th5_av3, inf)/norm(th3, inf);
rel_Linf_err_4_th = norm(th4-th5_av4, inf)/norm(th4, inf);

figure(1)
plot(log2(x), log10([rel_L2_err_1_th rel_L2_err_2_th rel_L2_err_3_th rel_L2_err_4_th]), 'kd-', 'Linewidth', 1.5, 'MarkerSize', 7.0)
hold on
plot(log2(x), log10([.5*rel_L2_err_1_th .5/4*rel_L2_err_1_th .5/16*rel_L2_err_1_th .5/64*rel_L2_err_1_th]), 'k-.', 'Linewidth', 1.5)
set(gca, 'Fontname', 'Times', 'fontsize', 16)
%title('Relative L^2 error');
xlabsize = get(gca, 'xlabel');
set(xlabsize,'string','log$_2$ N points','fontsize', 16, 'fontname', 'times')
ylabsize = get(gca, 'ylabel');
set(ylabsize,'string','log$_{10}$ L$^2$ error','fontsize', 16, 'fontname', 'times')

text(log2(.5*(x(1)+x(2))), double(log10(.7*rel_L2_err_1_th)), num2str(log2(rel_L2_err_1_th/rel_L2_err_2_th)));
text(log2(.5*(x(2)+x(3))), double(log10(.7*rel_L2_err_2_th)), num2str(log2(rel_L2_err_2_th/rel_L2_err_3_th)));
text(log2(.5*(x(3)+x(4))), double(log10(.7*rel_L2_err_3_th)), num2str(log2(rel_L2_err_3_th/rel_L2_err_4_th)));

fig=gcf;
fig.Color = 'white';
fig.InvertHardcopy = 'off';
filename = sprintf('../RKLM_Reference/Doc/paper_2019/figures/TravellingVortex/%s/%s_convergence_L2.eps', varstr, varstr);
print(filename, '-depsc')
export_fig(filename, '-eps')

figure(2)
plot(log2(x), log10([rel_Linf_err_1_th rel_Linf_err_2_th rel_Linf_err_3_th rel_Linf_err_4_th]), 'kd-', 'Linewidth', 1.5, 'MarkerSize', 7.0)
hold on
plot(log2(x), log10([.5*rel_Linf_err_1_th .5/4*rel_Linf_err_1_th .5/16*rel_Linf_err_1_th .5/64*rel_Linf_err_1_th]), 'k-.', 'Linewidth', 1.5)
set(gca, 'Fontname', 'Times', 'fontsize', 16)
%title('thity relative L^{\infty} error');  
xlabsize = get(gca, 'xlabel');
set(xlabsize,'string','log$_2$ N points','fontsize', 16, 'fontname', 'times')
ylabsize = get(gca, 'ylabel');
set(ylabsize,'string','log$_{10}$ L$^{\infty}$ error','fontsize', 16, 'fontname', 'times')

text(log2(.5*(x(1)+x(2))), double(log10(.7*rel_Linf_err_1_th)), num2str(log2(rel_Linf_err_1_th/rel_Linf_err_2_th)));
text(log2(.5*(x(2)+x(3))), double(log10(.7*rel_Linf_err_2_th)), num2str(log2(rel_Linf_err_2_th/rel_Linf_err_3_th)));
text(log2(.5*(x(3)+x(4))), double(log10(.7*rel_Linf_err_3_th)), num2str(log2(rel_Linf_err_3_th/rel_Linf_err_4_th)));

fig=gcf;
fig.Color = 'white';
fig.InvertHardcopy = 'off';
filename = sprintf('../RKLM_Reference/Doc/paper_2019/figures/TravellingVortex/%s/%s_convergence_Linfty.eps', varstr, varstr);
print(filename, '-depsc')
export_fig(filename, '-eps')
    