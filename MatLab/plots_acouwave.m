function plots_acouwave(varstr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% plots_acouwave(varstr)
% Produces plots from RKLM low Mach fluid dynamics code hdf output
% for the paper "A semi-implicit numerical model for small-to-planetary scale atmospheric
% dynamics", acoustic wave test case
%
%
% varstr    variable to plot
%           'p2_nodes', 'dp2_c', 'theta', 'dp2_nodes', 'dpdim', 'u',
%           'v', 'w', 'dY', 'rho', 'p', 'geopot', 'p2_c', 'rhoY'
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

title_true = 0;

kmin = 0;
kmax = 601;
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

folderstring = strcat('/home/tommaso/work/repos/RKLM_Reference/low_Mach_gravity_comp');

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

if strcmp(varstr,'p2_n')
    folderstr = 'p2_nodes';
    arraysize = [nnx nny];
end

set(0,'defaulttextinterpreter','latex')
scrsz = get(0,'ScreenSize');

figure1 = figure('Position',[1 2*scrsz(4)/3 scrsz(4)/2 1*scrsz(4)/3]);

figure2 = figure('Position',[scrsz(4)/2 2*scrsz(4)/3 scrsz(4)/2 1*scrsz(4)/3]);
title(strcat('horizontal slice at j = ',num2str(showslice_hor)))

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
        
        % Create filled contour
        figure(figure1)
        if strcmp(varstr, 'flux_rhou')
            contourf(x,z,th./th(:,1),no_of_lines,'LineColor',linecolor);
        elseif (strcmp(varstr, 'u') || strcmp(varstr, 'v') || strcmp(varstr, 'w'))
            contourf(x,z,th*velosc,no_of_lines,'LineColor',linecolor);
        else
            if diff_rel_to_bottom
                th = th-th(3,:);
            end
            contourf(x,z,th,no_of_lines,'LineColor',linecolor);
        end
        colormap viridis
        colorbar('FontSize',14,'FontName','Helvetica')
        if(k==1) % widen picture size
            pos=get(gca,'position');  % retrieve the current plot size value
            pos(3)=.95*pos(3);        % try increasing width and height 10%
            pos(4)=.95*pos(4);        % try increasing width and height 10%
            set(gca,'position',pos);  % write the new values
        end
        
        set(gca,'DataAspectRatio', aspect, 'FontSize',18,'FontName','Helvetica');
        axis tight;
        if title_true
            title(strcat(titlestr,kstr));
        end
        set(0,'defaulttextinterpreter','latex')
        
    end
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
    
    
    
    
    
    figure(figure2)
    plot(th(showslice_hor,:))
    hold on
    
     xlabel('x [m]','FontSize',18,'Interpreter','latex');
     ylabel('$\rho$ [$kg m^{-3}$]','FontSize',18,'Interpreter','latex');
     xlim([0 257])
     xticks([0 128 257])
     xticklabels({'0', '0.5', '1'})
     if varstr=='rho'
         ylim([0.97 1.03])
         yticks([0.98 1.0 1.02])
         xticklabels({'0.98', '1.0', '1.02'})
%     elseif varstr=='rhou'
%         ylim([-75 75])
%         yticks([-50 0 50])
     end
     grid on
    fig=gcf;
    fig.Color = 'white';
    fig.InvertHardcopy = 'off';
    filename = sprintf('../RKLM_Reference/Doc/paper_2019/figures/Acoustic-Wave/%s/%s_snapshot%d.eps', varstr, varstr, k);
    print(filename, '-depsc')
    export_fig(filename, '-eps')
    
    
end