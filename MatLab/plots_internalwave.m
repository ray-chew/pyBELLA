function plots_internalwave(varstr, ext)
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
%           'NH', 'H', 'P'
%
% Developed by R. Klein, FU Berlin, -2019
% Modified by T. Benacchio, Politecnico di Milano, 2019
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Hint:
% saving figures as .eps:     print(gcf, 'TestPlot', '-depsc');

addpath('/home/tommaso/work/code/matlab_packages/export_fig')
addpath('/home/tommaso/work/code/matlab_packages/Colormaps')

set(0,'DefaultFigureColor',[1 1 1])

extrafigno = 52;

showmode = 1;
separate_signs = 1;
linecolor      = 'default';  % 'k', 'default' ...
no_of_contours = 10;
no_of_lines    = 25;
show_increments = 0;
symmetry = 0;        % in {0,1}
showdummycells = 0;
showslice = 1;
diff_rel_to_bottom = 0;
print_eps = 1;

% th0 = -0.0015/300;
% dth = 5e-4/300;
% contour_values = [th0 th0+dth th0+2*dth th0+3*dth th0+4*dth th0+5*dth th0+6*dth th0+7*dth th0+8*dth th0+9*dth th0+10*dth];
dtheta = 0.5e-3/300;
contour_values = 10.0*[-5*dtheta, -4*dtheta, -3*dtheta, -2*dtheta, -dtheta, 0.0, dtheta, 2*dtheta, 3*dtheta, 4*dtheta, 5*dtheta];
%contour_values = [1.0001 1.0011 1.0022 1.0022 1.0033 1.0044 1.0055 1.0065];
%contour_values = linspace(-0.01,0.01,41) / 288.15;
title_true = 0;

kmin = 0;
kmax = 2;
dk   = 1;

if strcmp(ext, 'NH')
    scalefactor = 1.0;  % Skamarock-Klemp-1994 Fig.1
elseif strcmp(ext, 'H')
    scalefactor = 20.0;  % Skamarock-Klemp-1994 Fig.3
elseif strcmp(ext, 'P')
    scalefactor = 160.0;   % new, very long wave test
end

ncx =301;
ncy = 10;
L   = 300.0 * scalefactor;  %
x0  = 0.0;
H   = 10.0;  %
aspect = [L/H/3 1 1];
velosc = 100;  % velocity unit of RKLM code
showslice_hor = floor(ncy/2);
showslice_ver = floor(ncx/2);


% auxiliary adjustments of grid parameters
nnx = ncx+1;
dumsx = 2;
if ncy == 1
    nny = 1;
    dumsy = 0;
else
    nny = ncy+1;
    dumsy = 2;
end


rhoY_diff = 0;
rhoZ_diff = 0;
transp    = 0;

folderstring = strcat('/home/tommaso/work/repos/RKLM_Reference/hdf_output/InternalWave_', ext);

% for time series display
ts_name = strcat(folderstring, '/time_series.txt');
% [rho_ts,rhou_ts,rhov_ts,rhow_ts,rhoe_ts,rhoY_ts] = import_timeseries(ts_name, 2, nts);

% cell-centered fields
%%varstr = 'drhoY';  folderstr = 'drhoY'; titlestr = 'drhoY'; ndummy = 2; arraysize = [ncx ncy];
%varstr = 'T';  folderstr = 'T'; titlestr = 'T'; ndummy = 2; arraysize = [ncx ncy];
%varstr = 'dT';  folderstr = 'dT'; titlestr = 'dT'; ndummy = 2; arraysize = [ncx ncy];
%varstr = 'Y';  folderstr = 'Y'; titlestr = '\theta'; ndummy = 2; arraysize = [ncx ncy];
folderstr = varstr;
ndummy = 2;
arraysize = [ncx ncy];

%varstr = 'dY';  folderstr = 'dY'; titlestr = 'd\theta'; ndummy = 2; arraysize = [ncx ncy]; filledcontours = 0; fixed_contours = 1;

if strcmp(varstr,'dY')
    filledcontours = 0;
    fixed_contours = 1;
end

set(0,'defaulttextinterpreter','latex')
scrsz = get(0,'ScreenSize');

figure1 = figure('Position',[1 2*scrsz(4)/3 scrsz(4)/2 1*scrsz(4)/3]);

for k = kmin:dk:kmax
    kstr = num2str(k);
    if k < 10
        filestr = strcat(folderstring,'/',folderstr,'/',varstr,'_00',kstr,'.hdf');
    else
        if k < 100
            filestr = strcat(folderstring,'/',folderstr,'/',varstr,'_0',kstr,'.hdf');
        else
            filestr = strcat(folderstring,'/',folderstr,'/',varstr,'_',kstr,'.hdf');
        end
    end
    if nny == 1
        v = hdfread(filestr, '/Data-Set-2', 'Index', {1,1,arraysize(1)+dumsx*ndummy});
        if strcmp(varstr, 'flux_rhov')
            transpose(v);
        end
        velo = v;
        
        [ndum, nx] = size(velo);
        nx = nx - 4;
        dx = L/nx;
        
        x = linspace(x0 + 0.5*dx-(nx/2)*dx,x0 - 0.5*dx+(nx/2)*dx,nx);
        Yt = transpose(velo);
        th = Yt(3:1:nx+2);
        
        % 1D graphs
        plot(x,th, 'r+');
        
        % Create xlabel
        % xlabel('x [10 km]','FontSize',18,'FontName','Helvetica');
        xlabel('x [km]','FontSize',18,'FontName','Helvetica');
        
        % Create ylabel
        ylabel(titlestr,'FontSize',18,'FontName','Helvetica');
        
    else
        v = hdfread(filestr, '/Data-Set-2', 'Index', {[1  1],[1  1],[arraysize(1)+dumsx*ndummy  arraysize(2)+dumsy*ndummy]});
        
        if transp == 1
            v = transpose(v);
        end
        
        if rhoY_diff
            velo = v - rhoY_diff;
        else
            velo = v;
        end
        
        if showdummycells == 1
            [nx, nz] = size(velo);
            
            dx = L/(nx-4);
            dz = H/(nz-4);
            
            dth = 8.3333e-04;
            x = linspace(x0 + 0.5*dx-(nx/2)*dx,x0 - 0.5*dx+(nx/2)*dx,nx);
            z = linspace(0.5*dz,-0.5*dz+nz*dz,nz);
            Yt = transpose(velo);
            th = Yt;
        else
            [nx, nz] = size(velo);
            nx = nx - 4;
            nz = nz - 4;
            
            dx = L/nx;
            dz = H/nz;
            
            dth = 8.3333e-04;
            x = linspace(x0 + 0.5*dx-(nx/2)*dx,x0 - 0.5*dx+(nx/2)*dx,nx);
            z = linspace(0.5*dz,-0.5*dz+nz*dz,nz);
            Yt = transpose(velo);
            th = Yt(3:1:nz+2, 3:1:nx+2);
        end
        
        % Create unfilled contour
        figure(figure1)
        if fixed_contours
            if separate_signs == 1
                contour(x,z,max(0.0,th),contour_values,'LineColor','k','LineWidth',1.0);
                % contour(x,z,max(0.0,th),contour_values,'LineColor','k');
                hold
                contour(x,z,min(0.0,th),contour_values,'LineColor','k');
                % contour(x,z,min(0.0,th),contour_values,'LineColor',linecolor,'LineStyle','--');
                hold
            else
                contour(x,z,th,contour_values,'LineColor',linecolor);
            end
        elseif fixed_contour_step
            if separate_signs == 1
                contour(x,z,max(0.0,th), 'LevelStep', dtheta, 'LineColor','k');
                hold
                contour(x,z,min(0.0,th), 'LevelStep', dtheta, 'LineColor','k','LineStyle','--');
                hold
            else
                contour(x,z,th, 'LevelStep', dtheta, 'LineColor','k');
            end
        else
            if separate_signs == 1
                contour(x,z,max(0.0,th),no_of_contours,'LineColor','k');
                hold
                contour(x,z,min(0.0,th),no_of_contours,'LineColor','k','LineStyle','--');
                hold
            else
                contour(x,z,th,no_of_contours,'LineColor','k');
            end
        end
        
    end
    
    set(gca,'DataAspectRatio', aspect, 'FontSize',18,'FontName','Helvetica');
    axis tight;
    
    set(0,'defaulttextinterpreter','latex')
    set(0,'DefaultFigureColor',[1 1 1])
    xlabel('x [km]','FontSize',18,'Interpreter','latex');
    ylabel('z [km]','FontSize',18,'Interpreter','latex');
    ylim([0 10])
    yticks([0 2 4 6 8 10])
    
    switch ext
        case 'NH'
            xlim([-150 150])
            xticks([-150 -100 -50 0 50 100 150])
            xticklabels({'0', '50', '100', '150', '200', '250', '300'})
        case 'H'
            xlim([-3000 3000])
            xticks([-3000 -2000 -1000 0 1000 2000 3000])
            xticklabels({'0', '1000', '2000', '3000', '4000', '5000', '6000'})
        case 'P'
            xlabel('x [$10^3$ km]','FontSize',18,'Interpreter','latex');
            xlim([-24000 24000])
            xticks([-24000 -16000 -8000 0 8000 16000 24000])
            xticklabels({'0', '8', '16', '24', '32', '40', '48'})
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
    
    %        set(gca,'yticklabel',[], 'xticklabel', []) %Remove tick labels
    
    %         %%Get tick mark positions
    %         yTicks = get(gca,'ytick');
    %         xTicks = get(gca, 'xtick');
    %         ax = axis; %Get left most x-position
    %         HorizontalOffset = 0.1;
    %         %%Reset the ytick labels in desired font
    %         for i = 1:length(yTicks)
    %             %Create text box and set appropriate properties
    %             text(ax(1) - HorizontalOffset,yTicks(i),['$' num2str( yTicks(i)) '$'],...
    %                 'HorizontalAlignment','Right','interpreter', 'latex');
    %         end
    %         %%Reset the xtick labels in desired font
    %         minY = min(yTicks);
    %         verticalOffset = 0.2;
    %         for xx = 1:length(xTicks)
    %             %Create text box and set appropriate properties
    %             text(xTicks(xx), minY - verticalOffset, ['$' num2str( xTicks(xx)) '$'],...
    %                 'HorizontalAlignment','Right','interpreter', 'latex');
    %         end
    
    
end