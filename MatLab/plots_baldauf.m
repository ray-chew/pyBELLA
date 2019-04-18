function plots_baldauf(varstr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% plots_baldauf(varstr)
% Produces plots from RKLM low Mach fluid dynamics code hdf output
% for the paper "A semi-implicit numerical model for small-to-planetary scale atmospheric
% dynamics", large-scale internal wave test case with rotation
%
%
% varstr    variable to plot
%           'dT', 'u','v'
%
% Developed by R. Klein, FU Berlin, -2019
% Modified by T. Benacchio, Politecnico di Milano, 2019
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('./export_fig')

set(0,'DefaultFigureColor',[1 1 1])

kmin = 0;
kmax = 9;
dk   = 9;

scalefactor = 20.0;

ncx = 1201;
ncy = 80;
L   = 300.0 * scalefactor;  %
x0  = 0.0;
H   = 10.0;  %
aspect = [L/H/3 1 1];

if strcmp(varstr, 'dT')
    contour_values_p = linspace(0.0012, 0.006, 5);
    contour_values_m = linspace(-0.006, -0.0012, 5);
elseif strcmp(varstr, 'v')
    contour_values_p = linspace(0.0002, 0.0012, 6);
    contour_values_m = linspace(-0.0012, -0.0002, 6);
elseif strcmp(varstr, 'u')
    contour_values_p = linspace(0.002, 0.012, 6);
    contour_values_m = linspace(-0.012, -0.002, 6);
else
    disp('Error, please enter dT, u, or v as the variable to plot.')
    return;
end

% auxiliary adjustments of grid parameters
dumsx = 2;
dumsy = 2;

folderstring = strcat('../hdf_output/InternalWave_Baldauf_1200');

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
    th = Yt(3:1:nz+2, 3:1:nx+2);
    
    if strcmp(varstr, 'dT')
        dim_scale = 300.0;
    elseif strcmp(varstr, 'v')
        dim_scale = 10.0;
    elseif strcmp(varstr, 'u')
        dim_scale = 10.0;
    else
        disp('Error, please enter dT, u, or w as the variable to plot.')
        return;
    end
            
    % Create filled contour
    figure(figure1)
    [ccf1,hhf1]=contourf(x,z,dim_scale*th,[min(min(dim_scale*th)) contour_values_m contour_values_p max(max(dim_scale*th))], 'LineColor','k','LineWidth',1.0);
    set(hhf1,'LineColor','none')
    hold on
    contour(x,z,dim_scale*th,[contour_values_p max(max(dim_scale*th))],'LineColor','k','LineWidth',1.0);
    [cc1,h1]=contour(x,z,dim_scale*th,[min(min(dim_scale*th)) contour_values_m],'LineStyle', '--','LineColor','k','LineWidth',1.0);
    
    % Take all the info from the contourline output argument:
    i0 = 1;
    i2 = 1;
    while i0 <  length(cc1)
        i1 = i0+[1:cc1(2,i0)];
        zLevel(i2) = cc1(1,i0);
        hold on
        % And plot it with dashed lines:
        ph(i2) = plot(cc1(1,i1),cc1(2,i1),'k--','linewidth',1);
        i0 = i1(end)+1;
        i2 = i2+1;
    end
    % Scrap the contourlines:
    delete(h1)

    set(gca,'DataAspectRatio', aspect, 'FontSize',10,'FontName','Helvetica');
    axis tight;
    
    set(0,'defaulttextinterpreter','latex')
    set(0,'DefaultFigureColor',[1 1 1])
    xlim([-3000 3000])
    xticks([-3000 -2000 -1000 0 1000 2000 3000])
    xticklabels({'0', '1000', '2000', '3000', '4000', '5000', '6000'})
    xlabel('x [km]','FontSize',14,'Interpreter','latex');
    ylabel('z [km]','FontSize',14,'Interpreter','latex');
    yticks([0 2 4 6 8 10])
    colormap viridis
    colorbar('FontSize',10,'FontName','Helvetica');
        
    fig=gcf;
    fig.Color = 'white';
    fig.InvertHardcopy = 'off';
    filename = sprintf('../RKLM_Reference/Doc/paper_2019/figures/InternalWave_Baldauf_1200/%s/%s_snapshot%d.eps', varstr, varstr, k);
    print(filename, '-depsc')
    export_fig(filename, '-eps')
    hold off
end
