% Script collecting some convenient data views on the hdf4 output
% of the RKLM code family

% Hint: 
% saving figures as .eps:     print(gcf, 'TestPlot', '-depsc');

%clear all
%close all

extrafigno = 52;

%modelstr = '';
modelstr = 'comp';
%modelstr = 'psinc' ;  
%modelstr = 'psinc_w_adv_Ndt=3';
%modelstr = 'psinc_Ndt=3';
%modelstr = 'psinc_w_adv_Ndt=05';


%test_case = 'Baldaufs-Internal-Wave-Tests';
%test_case = 'Deep-Internal-Wave-Tests';
%test_case = 'Internal-Wave-Tests';
%test_case = 'Breaking-Wave-Tests';
%test_case = 'Rising-Bubble';
%test_case = 'Smolarkiewicz-Margolin-Breaking-Wave';
%test_case = 'Straka_50m';
test_case = 'Travelling-Vortex';
%test_case = 'Gresho-Vortex';
%test_case = 'Travelling-Hump';
%test_case = 'Acoustic-Wave';
%test_case = 'Schlutows-Wave';

showmode = 1;
separate_signs = 1;
filledcontours = 1;
linecolor      = 'default';  % 'k', 'default' ...
no_of_lines    = 25;
fixed_contours = 0;
fixed_contour_step = 0;
no_of_contours = 10;
show_increments = 0;
symmetry = 0;        % in {0,1}
symmetrytest = 0;
showdummycells = 0;
showslice = 1;
diff_rel_to_bottom = 0;
print_eps = 0;

% th0 = -0.0015/300;
% dth = 5e-4/300;
% contour_values = [th0 th0+dth th0+2*dth th0+3*dth th0+4*dth th0+5*dth th0+6*dth th0+7*dth th0+8*dth th0+9*dth th0+10*dth];
dtheta = 0.5e-3/300;
contour_values = 10.0*[-5*dtheta, -4*dtheta, -3*dtheta, -2*dtheta, -dtheta, 0.0, dtheta, 2*dtheta, 3*dtheta, 4*dtheta, 5*dtheta];
%contour_values = [1.0001 1.0011 1.0022 1.0022 1.0033 1.0044 1.0055 1.0065];
%contour_values = linspace(-0.01,0.01,41) / 288.15;
title_true = 1;

kmin = 0;
kmax = 601;
dk   = 1;

if strcmp(test_case, 'Baldaufs-Internal-Wave-Tests')
    scalefactor = 20.0;
    ncx = 151; % 301; 
    ncy =  10;  %  20;  
    L   = scalefactor*300.0;  % [km] 
    x0  = 0.0;
    H   = 10.0;               % [km]     
    aspect = [L/H/3 1 1];
    velosc = 10;  % velocity unit of RKLM code
    showslice_hor = floor(ncy/2);
    showslice_ver = floor(ncx/2);
elseif strcmp(test_case, 'Deep-Internal-Wave-Tests')
    ncx = 300; 
    ncy = 30;  
    L   = 40000.0;  % 
    x0  = 0.0;
    H   = 300.0;  %
    aspect = [L/H/3 1 1];
    velosc = 100;  % velocity unit of RKLM code
    showslice_hor = floor(ncy/2);
    showslice_ver = floor(ncx/2);
elseif strcmp(test_case, 'Internal-Wave-Tests')
    scalefactor = 1.0;  % Skamarock-Klemp-1994 Fig.1
    % scalefactor = 20.0;   % Skamarock-Klemp-1994 Fig.3
    %scalefactor = 160.0;   % new, very long wave test
    ncx = 301; 
    ncy = 10;  
    %ncx = 600; 
    %ncy = 80;  
    L   = 300.0 * scalefactor;  % 
    x0  = 0.0;
    H   = 10.0;  %
    aspect = [L/H/3 1 1];
    velosc = 100;  % velocity unit of RKLM code
    showslice_hor = floor(ncy/2);
    showslice_ver = floor(ncx/2);
elseif strcmp(test_case, 'Breaking-Wave-Tests')
    ncx = 240; 
    ncy = 120;  
    L   = 120;  % 
    x0  = 0.0;
    H   = 60.0;  %
    aspect = [L/H/3 1 1];
    velosc = 65.15;  % velocity unit of RKLM code
    showslice_hor = floor(ncy/2);
    showslice_ver = floor(ncx/2);
elseif strcmp(test_case, 'Acoustic-Wave')
    scalefactor = 1.0;
    ncx = 300; 
    ncy = 10;  
    L   = 300.0 * scalefactor;  % 
    x0  = 0.5*L;
    H   = 10.0;  %
    aspect = [16 1 1];
    velosc = 100;  % velocity unit of RKLM code
    showslice_hor = floor(ncy/2);
    showslice_ver = floor(ncx/2);
elseif strcmp(test_case, 'Rising-Bubble')
    ncx = 160;  
    ncy =  80;  
    L   = 2.0;  
    x0  = 1.0;
    H   = 1.0; 
    aspect = [1 1 1];
    velosc = 100;  % velocity unit of RKLM code
    showslice_hor = floor(ncy/2);
    showslice_ver = floor(ncx/2);
elseif strcmp(test_case, 'Travelling-Vortex')
    ncx = 48;  
    ncy = 48; 
    L   = 1.0;  
    x0  = 0.0;
    H   = 1.0; 
    aspect = [2 2 2];
    velosc = 1; 
    showslice_hor = floor(ncy/2);
    showslice_ver = floor(ncx/2);
elseif strcmp(test_case, 'Gresho-Vortex')
    ncx = 80;  
    ncy = 20; 
    L   = 4.0;  
    x0  = 0.5;
    H   = 1.0; 
    aspect = [2 2 2];
    velosc = 1; 
    showslice_hor = floor(ncy/2);
    showslice_ver = floor(ncx/2);
elseif strcmp(test_case, 'Travelling-Hump')
    ncx = 20;  
    ncy = 20; 
    L   = 1.0;  
    x0  = 0.5;
    H   = 1.0; 
    aspect = [2 2 2];
    velosc = 100;  % velocity unit of RKLM code
    showslice_hor = ncy/2;
    showslice_ver = floor(ncx/2);
elseif strcmp(test_case, 'Straka_50m')
    ncx = 1025;  
    ncy = 128;  
    L  = 51.2;  % 
    x0 = 0.0*L;
    H  = 6.4;  %
    aspect = [1 1 1];
    velosc = 100;  % velocity unit of RKLM code
    dtheta = 1.0/300.0;
    contour_values = linspace(-16.5*dtheta,-0.5*dtheta,17);
    showslice_hor = floor(ncy/3);
    showslice_ver = floor(ncx/2);
elseif strcmp(test_case, 'Smolarkiewicz-Margolin-Breaking-Wave')
    ncx = 240;  
    ncy = 120;  
    L   = 120.0;  % 
    x0  =  60.0;
    H   =  60.0;  %
    aspect = [1 1 1];
    velosc = 100;  % velocity unit of RKLM code
    showslice_hor = floor(ncy/2);
    showslice_ver = floor(ncx/2);
elseif strcmp(test_case, 'Schlutows-Wave')
    ncx = 120;  
    ncy = 120;  
    L   =  30.0;  % 
    x0  =   0.0;
    H   =  60.0;  %
    aspect = [1 1 1];
    velosc = 100;  % velocity unit of RKLM code
    showslice_hor = 5; % floor(ncy/2);
    showslice_ver = floor(ncx/2);
end
    
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


modelfigstr = strcat('  (',modelstr,')');

rhoY_diff = 0;
rhoZ_diff = 0;
transp    = 0;


folderstring = strcat('/Users/rupert/Documents/Computation/RKLM_Reference/low_Mach_gravity_',modelstr);
%folderstring = strcat('/home/tommaso/work/repos/RKLM_Reference/low_Mach_gravity_',modelstr);

% for time series display
ts_name = strcat(folderstring, '/time_series.txt');
% [rho_ts,rhou_ts,rhov_ts,rhow_ts,rhoe_ts,rhoY_ts] = import_timeseries(ts_name, 2, nts);

% cell-centered fields
%varstr = 'rho'; folderstr = 'rho'; titlestr = 'rho'; ndummy = 2; arraysize = [ncx ncy]; filledcontours = 1; fixed_contours = 0;
%varstr = 'p'; folderstr = 'p'; titlestr = 'p'; ndummy = 2; arraysize = [ncx ncy];
%varstr = 'S'; folderstr = 'S'; titlestr = 'S'; ndummy = 2; arraysize = [ncx ncy];
%varstr = 'rhoY';  folderstr = 'rhoY'; titlestr = 'rhoY'; ndummy = 2; arraysize = [ncx ncy]; rhoY_diff = 1;
%varstr = 'drhoY';  folderstr = 'drhoY'; titlestr = 'drhoY'; ndummy = 2; arraysize = [ncx ncy];
%varstr = 'T';  folderstr = 'T'; titlestr = 'T'; ndummy = 2; arraysize = [ncx ncy];
%varstr = 'dT';  folderstr = 'dT'; titlestr = 'dT'; ndummy = 2; arraysize = [ncx ncy]; filledcontours = 1; fixed_contours = 1;
%varstr = 'Y';  folderstr = 'Y'; titlestr = '\theta'; ndummy = 2; arraysize = [ncx ncy]; filledcontours = 0; no_of_contours = 50;
%varstr = 'dY';  folderstr = 'dY'; titlestr = 'd\theta'; ndummy = 2; arraysize = [ncx ncy]; filledcontours = 1; fixed_contours = 1;
%varstr = 'buoy';  folderstr = 'buoy'; titlestr = 'buoy'; ndummy = 2; arraysize = [ncx ncy];
%varstr = 'rhoZp';  folderstr = 'rhoZp'; titlestr = 'rhoZp'; ndummy = 2; arraysize = [ncx ncy];
%varstr = 'rhoZB';  folderstr = 'rhoZB'; titlestr = 'rhoZB'; ndummy = 2; arraysize = [ncx ncy];
%varstr = 'rhoZ';  folderstr = 'rhoZ'; titlestr = 'rhoZ'; ndummy = 2; arraysize = [ncx ncy];
%varstr = 'rhou';  folderstr = 'rhou'; titlestr = 'rhou'; ndummy = 2; arraysize = [ncx ncy];
%varstr = 'rhov';  folderstr = 'rhov'; titlestr = 'rhov'; ndummy = 2; arraysize = [ncx ncy];
%varstr = 'rhow';  folderstr = 'rhow'; titlestr = 'rhow'; ndummy = 2; arraysize = [ncx ncy];
%varstr = 'u';  folderstr = 'u'; titlestr = 'u'; ndummy = 2; arraysize = [ncx ncy]; symmetry = -1*symmetry;
%varstr = 'v';  folderstr = 'v'; titlestr = 'v'; ndummy = 2; arraysize = [ncx ncy]; symmetry = -1*symmetry; filledcontours = 1;
%varstr = 'w';  folderstr = 'w'; titlestr = 'w'; ndummy = 2; arraysize = [ncx ncy];
%varstr = 'vortz';  folderstr = 'vortz'; titlestr = 'vortz'; ndummy = 2; arraysize = [nnx nny]; filledcontours = 1; fixed_contours = 0;
%varstr = 'qv';  folderstr = 'qv'; titlestr = 'qv'; ndummy = 2; arraysize = [ncx ncy];
%varstr = 'qc';  folderstr = 'qc'; titlestr = 'qc'; ndummy = 2; arraysize = [ncx ncy];
%varstr = 'qr';  folderstr = 'qr'; titlestr = 'qr'; ndummy = 2; arraysize = [ncx ncy];
%varstr = 'flux_rhou';  folderstr = 'fluxes'; titlestr = 'flux_rhou'; ndummy = 2; arraysize = [ncx+1 ncy]; symmetry = -1*symmetry;
%varstr = 'flux_rhov';  folderstr = 'fluxes'; titlestr = 'flux_rhov'; ndummy = 2; arraysize = [ncy+1 ncx]; symmetry = -1*symmetry;

%varstr = 'p2_c';  folderstr = 'p2_c'; titlestr = '\pi'; ndummy = 2; arraysize = [ncx ncy];
%varstr = 'dp2_c';  folderstr = 'dp2_c'; titlestr = 'd\pi'; ndummy = 2; arraysize = [ncx ncy];
%varstr = 'dpdim';  folderstr = 'dpdime'; titlestr = 'dp [Pa]'; ndummy = 2; arraysize = [ncx ncy];
%varstr = 'rhs_cells';  folderstr = 'rhs_cells'; titlestr = 'rhs_c'; ndummy = 2; arraysize = [ncx ncy];
%varstr = 'rhs_cells_prec';  folderstr = 'rhs_cells'; titlestr = 'rhs_c_prec'; ndummy = 2; arraysize = [ncx ncy];
%varstr = 'lap_cells';  folderstr = 'lap_cells'; titlestr = 'lap_c';    ndummy = 2; arraysize = [ncx ncy];
%varstr = 'dP_cells';  folderstr = 'rhs_cells'; titlestr = 'dP_c'; ndummy = 2; arraysize = [ncx ncy];

varstr = 'p2_n';  folderstr = 'p2_nodes'; titlestr = '\pi_n';    ndummy = 2; arraysize = [nnx nny];
%varstr = 'dp2_n';  folderstr = 'dp2_nodes'; titlestr = 'd\pi_n';    ndummy = 2; arraysize = [nnx nny];
%varstr = 'rhs_nodes';  folderstr = 'rhs_nodes'; titlestr = 'rhs_n';    ndummy = 2; arraysize = [nnx nny];
%varstr = 'rhs_nodes_prec';  folderstr = 'rhs_nodes'; titlestr = 'rhs_n_prec';    ndummy = 2; arraysize = [nnx nny];
%varstr = 'lap_nodes';  folderstr = 'lap_nodes'; titlestr = 'lap_n';    ndummy = 2; arraysize = [nnx nny];

%varstr = 'flux_x';  folderstr = 'flux_x'; titlestr = 'flux_x'; ndummy = 2; arraysize = [ncx+1 ncy]; symmetry = -1*symmetry;
%varstr = 'flux_y';  folderstr = 'flux_y'; titlestr = 'flux_y'; ndummy = 2; arraysize = [ncy+1 ncx]; symmetry = -1*symmetry; transp = 1;


scrsz = get(0,'ScreenSize');
figure1 = figure('Position',[1 2*scrsz(4)/3 scrsz(4)/2 1*scrsz(4)/3]);
if abs(symmetry) == 1 
    figure2 = figure('Position',[scrsz(4)/2 0 scrsz(3)/1 scrsz(4)/2.5]);
end

if title_true == 0
    texthandle = annotation('textbox', [0.45 0.95 0.2 0.05]);
    set(texthandle, 'String', strcat(titlestr,modelfigstr,' 2D'), 'FontSize', 18, 'FontName', 'Helvetica', 'LineStyle', 'none');
end
%set(texthandle, 'String', titlestr, 'FontSize', 14, 'FontName', 'Optima');

if showslice_hor
    figure3 = figure('Position',[scrsz(4)/2 2*scrsz(4)/3 scrsz(4)/2 1*scrsz(4)/3]);
    title(strcat('horizontal slice at j = ',num2str(showslice_hor)))
    set(gca,'FontSize',18,'FontName','Helvetica');
end
if showslice_ver
    figure4 = figure('Position',[2*scrsz(4)/2 2*scrsz(4)/3 scrsz(4)/2 1*scrsz(4)/3]);
    title(strcat('vertical slice at i = ',num2str(showslice_ver)))
    set(gca,'FontSize',18,'FontName','Helvetica');
end

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
        
        % Create contour
        if filledcontours
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
                disp(min(min(th))*300)
            end
            colormap Jet
            colorbar('FontSize',14,'FontName','Helvetica')
            if(k==1) % widen picture size
                pos=get(gca,'position');  % retrieve the current plot size value
                pos(3)=.95*pos(3);        % try increasing width and height 10% 
                pos(4)=.95*pos(4);        % try increasing width and height 10% 
                set(gca,'position',pos);  % write the new values
            end
            if print_eps
                filename = sprintf('./results/%s_evol/%s_snapshot%d.eps', varstr, varstr, k);
                print(filename, '-depsc')            
            end
        else
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
            if print_eps
                filename = sprintf('./results/%s_evol/%s_snapshot%d.eps', varstr, varstr, k);
                print(filename, '-depsc')
            end
        end
                
        set(gca,'DataAspectRatio', aspect, 'FontSize',18,'FontName','Helvetica');
        axis tight;
        if title_true
            title(strcat(titlestr,kstr));
        end
        % Create xlabel
        % xlabel('x [10 km]','FontSize',18,'FontName','Helvetica');
        xlabel('x [km]','FontSize',18,'FontName','Helvetica');
        
        % Create ylabel
        ylabel('z [km]','FontSize',18,'FontName','Helvetica');

        
        
        
        if show_increments
            if k>kmin
                figure(80)
                contourf(x,z,th-th_old,no_of_lines,'LineColor',linecolor);
                colormap Jet
                colorbar('FontSize',14,'FontName','Helvetica')
            end
            th_old = th;
        end

        if abs(symmetry) == 1
            figure(figure2);
            NonSymmAmpl = CheckSymmetry(th,x,z,k,symmetry);
            colormap Jet
            colorbar('FontSize',14,'FontName','Helvetica')
            %NonSymmAmpl;
        end

    end
         
    if symmetrytest == 1
        SymmetryTests(transpose(th),55,[4 1 1]);
        figure(figure1)
    end
    
    if (strcmp(varstr, 'u') || strcmp(varstr, 'v') || strcmp(varstr, 'w'))
        th = th*velosc;
    end
    if showslice_hor
        figure(figure3)
        hold
        plot(th(showslice_hor,:))
        if print_eps
            filename = sprintf('./results/%s_evol/%s_snapshot%d_cut_hor.eps', varstr, varstr, k);
            print(filename, '-depsc')
        else
            hold
        end
    end
    if showslice_ver
        figure(figure4)
        hold    
        plot(th(:,showslice_hor))
        if print_eps
            filename = sprintf('./results/%s_evol/%s_snapshot%d_cut_ver.eps', varstr, varstr, k);
            print(filename, '-depsc')
        else
            hold
        end
    end
    
    pause
end


