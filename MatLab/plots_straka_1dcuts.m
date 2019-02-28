function plots_straka_1dcuts(varstr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% plots_straka_1dcuts(varstr, ext)
% Produces plots from RKLM low Mach fluid dynamics code hdf output
% for the paper "A semi-implicit numerical model for small-to-planetary scale atmospheric
% dynamics",density current test case
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

addpath('/home/tommaso/work/code/matlab_packages/export_fig')
addpath('/home/tommaso/work/code/matlab_packages/Colormaps')

set(0,'DefaultFigureColor',[1 1 1])

linecolor      = 'default';  % 'k', 'default' ...

kmin = 0;
kmax = 3;
dk   = 1;


L  = 51.2;  
x0 = 0.0*L;
H  = 6.4;  
aspect = [1 1 1];
dtheta = 1.0/300.0;
contour_values = linspace(-16.5*dtheta,-0.5*dtheta,16);

% auxiliary adjustments of grid parameters
dumsx = 2;
dumsy = 2;

folderstring1 = strcat('/home/tommaso/work/repos/RKLM_Reference/hdf_output/Straka_3D_400m');
folderstring2 = strcat('/home/tommaso/work/repos/RKLM_Reference/hdf_output/Straka_3D_200m');
folderstring3 = strcat('/home/tommaso/work/repos/RKLM_Reference/hdf_output/Straka_3D_100m');
folderstring4 = strcat('/home/tommaso/work/repos/RKLM_Reference/hdf_output/Straka_3D_50m');

% cell-centered fields
folderstr = varstr;
ndummy = 2;

set(0,'defaulttextinterpreter','latex')
scrsz = get(0,'ScreenSize');

figure1 = figure('Position',[1 2*scrsz(4)/3 scrsz(4)/2 1*scrsz(4)/3]);

k = kmax; % Just plotting at end time
kstr = num2str(k);
filestr1 = strcat(folderstring1,'/',folderstr,'/',varstr,'_00',kstr,'.hdf');
filestr2 = strcat(folderstring2,'/',folderstr,'/',varstr,'_00',kstr,'.hdf');
filestr3 = strcat(folderstring3,'/',folderstr,'/',varstr,'_00',kstr,'.hdf');
filestr4 = strcat(folderstring4,'/',folderstr,'/',varstr,'_00',kstr,'.hdf');

ncx = 129;
ncy = 16;
arraysize = [ncx ncy];
v1 = hdfread(filestr1, '/Data-Set-2', 'Index', {[1  1],[1  1],[arraysize(1)+dumsx*ndummy  arraysize(2)+dumsy*ndummy]});
ncx = 257;
ncy = 32;
arraysize = [ncx ncy];
v2 = hdfread(filestr2, '/Data-Set-2', 'Index', {[1  1],[1  1],[arraysize(1)+dumsx*ndummy  arraysize(2)+dumsy*ndummy]});
ncx = 513;
ncy = 64;
arraysize = [ncx ncy];
v3 = hdfread(filestr3, '/Data-Set-2', 'Index', {[1  1],[1  1],[arraysize(1)+dumsx*ndummy  arraysize(2)+dumsy*ndummy]});
ncx = 1025;
ncy = 128;
arraysize = [ncx ncy];
v4 = hdfread(filestr4, '/Data-Set-2', 'Index', {[1  1],[1  1],[arraysize(1)+dumsx*ndummy  arraysize(2)+dumsy*ndummy]});

[nx1, nz1] = size(v1);
nx1 = nx1 - 4;
nz1 = nz1 - 4;

dx1 = L/nx1;
dz1 = H/nz1;

x1 = linspace(x0 + 0.5*dx1-(nx1/2)*dx1,x0 - 0.5*dx1+(nx1/2)*dx1,nx1);
z1 = linspace(0.5*dz1,-0.5*dz1+nz1*dz1,nz1);
Yt1 = transpose(v1);
th1 = Yt1(3:1:nz1+2, 3:1:nx1+2);


[nx2, nz2] = size(v2);
nx2 = nx2 - 4;
nz2 = nz2 - 4;

dx2 = L/nx2;
dz2 = H/nz2;

x2 = linspace(x0 + 0.5*dx2-(nx2/2)*dx2,x0 - 0.5*dx2+(nx2/2)*dx2,nx2);
z2 = linspace(0.5*dz2,-0.5*dz2+nz2*dz2,nz2);
Yt2 = transpose(v2);
th2 = Yt2(3:1:nz2+2, 3:1:nx2+2);

[nx3, nz3] = size(v3);
nx3 = nx3 - 4;
nz3 = nz3 - 4;

dx3 = L/nx3;
dz3 = H/nz3;

x3 = linspace(x0 + 0.5*dx3-(nx3/2)*dx3,x0 - 0.5*dx3+(nx3/2)*dx3,nx3);
z3 = linspace(0.5*dz3,-0.5*dz3+nz3*dz3,nz3);
Yt3 = transpose(v3);
th3 = Yt3(3:1:nz3+2, 3:1:nx3+2);

[nx4, nz4] = size(v4);
nx4 = nx4 - 4;
nz4 = nz4 - 4;

dx4 = L/nx4;
dz4 = H/nz4;

x4 = linspace(x0 + 0.5*dx4-(nx4/2)*dx4,x0 - 0.5*dx4+(nx4/2)*dx4,nx4);
z4 = linspace(0.5*dz4,-0.5*dz4+nz4*dz4,nz4);
Yt4 = transpose(v4);
th4 = Yt4(3:1:nz4+2, 3:1:nx4+2);

horx_slice = 3/16; %z=1200/6400 m
hor_xlim=40.6/51.2; % xlim[0 15 km]

% Create filled contour
figure(1)
plot(x1(floor(end/2+1):(floor(hor_xlim*end)+1)), ...
    0.5*(th1(floor(horx_slice*end), floor(end/2+1):(floor(hor_xlim*end)+1))+...
    th1(floor(horx_slice*end+1), floor(end/2+1):(floor(hor_xlim*end)+1)))*300, ...
    'k--', 'LineWidth',1.5,  'DisplayName', '400 m');
hold on
%figure(2)
plot(x2(floor(end/2+1):(floor(hor_xlim*end)+1)),...
    0.5*(th2(floor(horx_slice*end), floor(end/2+1):(floor(hor_xlim*end)+1))+...
    th2(floor(horx_slice*end+1), floor(end/2+1):(floor(hor_xlim*end)+1)))*300,...
    'r-.', 'LineWidth',1.5,  'DisplayName', '200 m');
%figure(3)
plot(x3(floor(end/2+1):(floor(hor_xlim*end)+1)),...
    0.5*(th3(floor(horx_slice*end), floor(end/2+1):(floor(hor_xlim*end)+1))+...
    th3(floor(horx_slice*end+1), floor(end/2+1):(floor(hor_xlim*end)+1)))*300,...
    'bo', 'LineWidth',1.5,  'DisplayName', '100 m');
%figure(4)
plot(x4(floor(end/2+1):(floor(hor_xlim*end)+1)), ...
    0.5*(th4(floor(horx_slice*end), floor(end/2+1):(floor(hor_xlim*end)+1))...
    +th4(floor(horx_slice*end+1), floor(end/2+1):(floor(hor_xlim*end)+1)))*300, ...
    'm-', 'LineWidth',1.5,  'DisplayName', '50 m');

set(gca,'DataAspectRatio', aspect, 'FontSize',18,'FontName','Helvetica');
axis tight;

set(0,'defaulttextinterpreter','latex')
set(0,'DefaultFigureColor',[1 1 1])
xlim([0 15])
xticks([0 5 10 15])
xlabel('x [km]','FontSize',18,'Interpreter','latex');
ylabel('$\theta'' $ [km]','FontSize',18,'Interpreter','latex');
%ylim([-5 1])
%yticks([-5 -4 -3 -2 -1 0 1])
grid on
legend('Location', 'SouthWest')

fig=gcf;
fig.Color = 'white';
fig.InvertHardcopy = 'off';
filename = sprintf('../RKLM_Reference/Doc/paper_2019/figures/Straka/%s/%s_snapshot1dcut.eps', varstr, varstr);
print(filename, '-depsc')
export_fig(filename, '-eps')
hold off
   
