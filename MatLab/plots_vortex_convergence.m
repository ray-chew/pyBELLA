function plots_vortex_convergence(varstr, np, conv_type)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% plots_vortex_convergence(varstr)
% Produces convergence plots from RKLM low Mach fluid dynamics code hdf output
% for the paper "A semi-implicit numerical model for small-to-planetary scale atmospheric
% dynamics", rotating Gresho vortex case. Assumes fine grid as exact
% solution
%
%
% varstr     variable to plot
%            'p2_nodes', 'dp2_c', 'theta', 'dp2_nodes', 'dpdim', 'u',
%            'v', 'w', 'dY', 'rho', 'p', 'geopot', 'p2_c', 'rhoY'
%
% np         number of points from which the grid is refined
%            e.g. 48 will check the errors on 48x48, 96x96, etc...
%
% conv_type  convergence kind
%            'self' takes fine-grid solution as exact and interpolates it
%            to coarser grids to compare it with finer grid solutions
%            'wrt_t0' takes initial solution as reference
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
kmax = 1;

%x = [48 96 192 384];
%x = [64 128 256 512];


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


% npxnp

ncx = np;
ncy = np;
nnx = ncy+1;
nny = ncy+1;

if strcmp(varstr, 'p2_n')
    arraysize = [nnx nny];
else 
    arraysize = [ncx ncy];
end

folderstring = strcat('../hdf_output/TravellingVortex_3D_64');
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
if strcmp(varstr, 'p2_n') && strcmp(conv_type, 'self')
    th1=th1(2:nz, 2:nx+2);
end

if strcmp(conv_type, 'wrt_t0')
    % Loading initial solution if convergence study with respect to initial
    % data
    filestr = strcat(folderstring,'/',folderstr,'/',varstr,'_000.hdf');
    v = hdfread(filestr, '/Data-Set-2', 'Index', {[1  1],[1  1],[arraysize(1)+dumsx*ndummy  arraysize(2)+dumsy*ndummy]});
    Yt = transpose(v);
    % Copy-pasting the first and last rows to get a square matrix
    th1_0(:, 1)      = Yt(3:1:nz+2, 3);
    th1_0(:, 2:nx+1) = Yt(3:1:nz+2, 3:1:nx+2);
    th1_0(:, nx+2)   = Yt(3:1:nz+2, nx+2);
    
end

% first refinement
ncx = 2*ncx;
ncy = 2*ncy;
nnx = ncy+1;
nny = ncy+1;

if strcmp(varstr, 'p2_n')
    arraysize = [nnx nny];
else
    arraysize = [ncx ncy];
end


folderstring = strcat('../hdf_output/TravellingVortex_3D_128');
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
if strcmp(varstr, 'p2_n') && strcmp(conv_type, 'self')
    th2=th2(2:nz, 2:nx+2);
end


if strcmp(conv_type, 'wrt_t0')
    % Loading initial solution if convergence study with respect to initial
    % data
    filestr = strcat(folderstring,'/',folderstr,'/',varstr,'_000.hdf');
    v = hdfread(filestr, '/Data-Set-2', 'Index', {[1  1],[1  1],[arraysize(1)+dumsx*ndummy  arraysize(2)+dumsy*ndummy]});
    Yt = transpose(v);
    % Copy-pasting the first and last rows to get a square matrix
    th2_0(:, 1)      = Yt(3:1:nz+2, 3);
    th2_0(:, 2:nx+1) = Yt(3:1:nz+2, 3:1:nx+2);
    th2_0(:, nx+2)   = Yt(3:1:nz+2, nx+2);
    
end




% Second refinement
ncx = 2*ncx;
ncy = 2*ncy;
nnx = ncy+1;
nny = ncy+1;

if strcmp(varstr, 'p2_n')
    arraysize = [nnx nny];
else
    arraysize = [ncx ncy];
end
folderstring = strcat('../hdf_output/TravellingVortex_3D_256');
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
if strcmp(varstr, 'p2_n') && strcmp(conv_type, 'self')
    th3=th3(2:nz, 2:nx+2);
end


if strcmp(conv_type, 'wrt_t0')
    % Loading initial solution if convergence study with respect to initial
    % data
    filestr = strcat(folderstring,'/',folderstr,'/',varstr,'_000.hdf');
    v = hdfread(filestr, '/Data-Set-2', 'Index', {[1  1],[1  1],[arraysize(1)+dumsx*ndummy  arraysize(2)+dumsy*ndummy]});
    Yt = transpose(v);
    % Copy-pasting the first and last rows to get a square matrix
    th3_0(:, 1)      = Yt(3:1:nz+2, 3);
    th3_0(:, 2:nx+1) = Yt(3:1:nz+2, 3:1:nx+2);
    th3_0(:, nx+2)   = Yt(3:1:nz+2, nx+2);
    
end



% third refinement
ncx = 2*ncx;
ncy = 2*ncy;
nnx = ncy+1;
nny = ncy+1;

if strcmp(varstr, 'p2_n')
    arraysize = [nnx nny];
else
    arraysize = [ncx ncy];
end

folderstring = strcat('../hdf_output/TravellingVortex_3D_512');
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
if strcmp(varstr, 'p2_n') && strcmp(conv_type, 'self')
    th4=th4(2:nz, 2:nx+2);
end

if strcmp(conv_type, 'wrt_t0')
    % Loading initial solution if convergence study with respect to initial
    % data
    filestr = strcat(folderstring,'/',folderstr,'/',varstr,'_000.hdf');
    v = hdfread(filestr, '/Data-Set-2', 'Index', {[1  1],[1  1],[arraysize(1)+dumsx*ndummy  arraysize(2)+dumsy*ndummy]});
    Yt = transpose(v);
    % Copy-pasting the first and last rows to get a square matrix
    th4_0(:, 1)      = Yt(3:1:nz+2, 3);
    th4_0(:, 2:nx+1) = Yt(3:1:nz+2, 3:1:nx+2);
    th4_0(:, nx+2)   = Yt(3:1:nz+2, nx+2);
    
end



% fourth refinement (this is the exact solution if conv_type='self')
ncx = 2*ncx;
ncy = 2*ncy;
nnx = ncy+1;
nny = ncy+1;

if strcmp(varstr, 'p2_n')
    arraysize = [nnx nny];
else
    arraysize = [ncx ncy];
end

folderstring = strcat('../hdf_output/TravellingVortex_3D_1024');
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
if strcmp(varstr, 'p2_n') && strcmp(conv_type, 'self')
    th5=th5(2:nz, 2:nx+2);
end

set(0,'defaulttextinterpreter','latex')
scrsz = get(0,'ScreenSize');



if strcmp(conv_type, 'wrt_t0')
    % Loading initial solution if convergence study with respect to initial
    % data
    filestr = strcat(folderstring,'/',folderstr,'/',varstr,'_000.hdf');
    v = hdfread(filestr, '/Data-Set-2', 'Index', {[1  1],[1  1],[arraysize(1)+dumsx*ndummy  arraysize(2)+dumsy*ndummy]});
    Yt = transpose(v);
    % Copy-pasting the first and last rows to get a square matrix
    th5_0(:, 1)      = Yt(3:1:nz+2, 3);
    th5_0(:, 2:nx+1) = Yt(3:1:nz+2, 3:1:nx+2);
    th5_0(:, nx+2)   = Yt(3:1:nz+2, nx+2);
    
end


if strcmp(conv_type, 'self')
    % ========================================================
    % Creating coarse-grid versions of fine-grid solution
    % ========================================================
    
    
    for j = 1:2:size(th5, 2)-1
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
    
    x_vec = [np 2*np 4*np 8*np];
    
    err_vec_2   = [rel_L2_err_1_th rel_L2_err_2_th rel_L2_err_3_th rel_L2_err_4_th];
    err_vec_inf = [rel_Linf_err_1_th rel_Linf_err_2_th rel_Linf_err_3_th rel_Linf_err_4_th];    

    scnd_ord = [.5*rel_L2_err_1_th .5/4*rel_L2_err_1_th .5/16*rel_L2_err_1_th .5/64*rel_L2_err_1_th];

else
    
    rel_L2_err_1_th = norm(th1-th1_0, 2)/norm(th1, 2);
    rel_L2_err_2_th = norm(th2-th2_0, 2)/norm(th2, 2);
    rel_L2_err_3_th = norm(th3-th3_0, 2)/norm(th3, 2);
    rel_L2_err_4_th = norm(th4-th4_0, 2)/norm(th4, 2);
    rel_L2_err_5_th = norm(th5-th5_0, 2)/norm(th5, 2);
    
    rel_Linf_err_1_th = norm(th1-th1_0, inf)/norm(th1, inf);
    rel_Linf_err_2_th = norm(th2-th2_0, inf)/norm(th2, inf);
    rel_Linf_err_3_th = norm(th3-th3_0, inf)/norm(th3, inf);
    rel_Linf_err_4_th = norm(th4-th4_0, inf)/norm(th4, inf);
    rel_Linf_err_5_th = norm(th5-th5_0, inf)/norm(th5, inf);

    x_vec = [np 2*np 4*np 8*np 16*np];
    
    
    
    err_vec_2   = [rel_L2_err_1_th rel_L2_err_2_th rel_L2_err_3_th rel_L2_err_4_th rel_L2_err_5_th];
    err_vec_inf = [rel_Linf_err_1_th rel_Linf_err_2_th rel_Linf_err_3_th rel_Linf_err_4_th rel_Linf_err_5_th];    

    scnd_ord = [.5*rel_L2_err_1_th .5/4*rel_L2_err_1_th .5/16*rel_L2_err_1_th .5/64*rel_L2_err_1_th .5/256*rel_L2_err_1_th];

    
end

figure1 = figure('Position',[1 2*scrsz(4)/3 scrsz(4)/2 1*scrsz(4)/3]);

figure(figure1)
plot(log2(x_vec), log10(err_vec_2), 'kd-', 'Linewidth', 1.5, 'MarkerSize', 7.0)
hold on
plot(log2(x_vec), log10(scnd_ord), 'k-.', 'Linewidth', 1.5)
set(gca, 'Fontname', 'Times', 'fontsize', 16)
%title('Relative L^2 error');
xlabsize = get(gca, 'xlabel');
set(xlabsize,'string','log$_2$ N points','fontsize', 16, 'fontname', 'times')
ylabsize = get(gca, 'ylabel');
set(ylabsize,'string','log$_{10}$ L$^2$ error','fontsize', 16, 'fontname', 'times')

text(log2(.5*(x_vec(1)+x_vec(2))), double(log10(.7*rel_L2_err_1_th)), num2str(log2(rel_L2_err_1_th/rel_L2_err_2_th)), 'fontsize', 16);
text(log2(.5*(x_vec(2)+x_vec(3))), double(log10(.7*rel_L2_err_2_th)), num2str(log2(rel_L2_err_2_th/rel_L2_err_3_th)), 'fontsize', 16);
text(log2(.5*(x_vec(3)+x_vec(4))), double(log10(.7*rel_L2_err_3_th)), num2str(log2(rel_L2_err_3_th/rel_L2_err_4_th)), 'fontsize', 16);
if strcmp(conv_type, 'wrt_t0')
   text(log2(.5*(x_vec(4)+x_vec(5))), double(log10(.7*rel_L2_err_4_th)), num2str(log2(rel_L2_err_4_th/rel_L2_err_5_th)), 'fontsize', 16);
end
axis square

fig=gcf;
fig.Color = 'white';
fig.InvertHardcopy = 'off';
filename = sprintf('../RKLM_Reference/Doc/paper_2019/figures/TravellingVortex/%s/%s_convergence_L2_%i_%s.eps', varstr, varstr, np, conv_type);
print(filename, '-depsc')
export_fig(filename, '-eps')
hold off
close;

figure2 = figure('Position',[1 2*scrsz(4)/3 scrsz(4)/2 1*scrsz(4)/3]);

figure(figure2)
plot(log2(x_vec), log10(err_vec_inf), 'kd-', 'Linewidth', 1.5, 'MarkerSize', 7.0)
hold on
plot(log2(x_vec), log10(scnd_ord), 'k-.', 'Linewidth', 1.5)
set(gca, 'Fontname', 'Times', 'fontsize', 16)
%title('thity relative L^{\infty} error');
xlabsize = get(gca, 'xlabel');
set(xlabsize,'string','log$_2$ N points','fontsize', 16, 'fontname', 'times')
ylabsize = get(gca, 'ylabel');
set(ylabsize,'string','log$_{10}$ L$^{\infty}$ error','fontsize', 16, 'fontname', 'times')

text(log2(.5*(x_vec(1)+x_vec(2))), double(log10(.7*rel_Linf_err_1_th)), num2str(log2(rel_Linf_err_1_th/rel_Linf_err_2_th)), 'fontsize', 16);
text(log2(.5*(x_vec(2)+x_vec(3))), double(log10(.7*rel_Linf_err_2_th)), num2str(log2(rel_Linf_err_2_th/rel_Linf_err_3_th)), 'fontsize', 16);
text(log2(.5*(x_vec(3)+x_vec(4))), double(log10(.7*rel_Linf_err_3_th)), num2str(log2(rel_Linf_err_3_th/rel_Linf_err_4_th)), 'fontsize', 16);
if strcmp(conv_type, 'wrt_t0')
   text(log2(.5*(x_vec(4)+x_vec(5))), double(log10(.7*rel_Linf_err_4_th)), num2str(log2(rel_Linf_err_4_th/rel_Linf_err_5_th)), 'fontsize', 16);
end
axis square

fig=gcf;
fig.Color = 'white';
fig.InvertHardcopy = 'off';
filename = sprintf('../RKLM_Reference/Doc/paper_2019/figures/TravellingVortex/%s/%s_convergence_Linfty_%i_%s.eps', varstr, varstr, np, conv_type);
print(filename, '-depsc')
export_fig(filename, '-eps')
