
% Hint: 
% saving figures as .eps:     print(gcf, 'TestPlot', '-depsc');

%modelstr = '';
%modelstr = 'comp';
modelstr = 'psinc';
%modelstr = 'psinc_w_adv_Ndt=3';
%modelstr = 'psinc_Ndt=3';
%modelstr = 'psinc_w_adv_Ndt=05';

%test_case = 'Equatorial-Long-Wave';
%test_case = 'Internal-Wave-Long-Wave';
%test_case = 'Internal-Wave-Strong-Strat';
%test_case = 'Skamarock-Klemp-Internal-Wave';
%test_case = 'Rising-Bubble';
%test_case = 'Smolarkiewicz-Margolin-Breaking-Wave';
%test_case = 'Straka';
test_case = 'Travelling-Vortex';
%test_case = 'Advection';

slice = 'full3D'; % options:  'xy' 'zy' 'full3D'
showmode = 1;
separate_signs = 1;
filledcontours = 1;
fixed_contours = 1;
fixed_contour_step = 0;
no_of_contours = 10;
show_increments = 0;
symmetry = 0;        % in {0,1}
symmetrytest = 0;
showdummycells = 0;

% th0 = -0.0015/300;
% dth = 5e-4/300;
% contour_values = [th0 th0+dth th0+2*dth th0+3*dth th0+4*dth th0+5*dth th0+6*dth th0+7*dth th0+8*dth th0+9*dth th0+10*dth];
dtheta = 0.5e-3/300;
contour_values = [-5*dtheta, -4*dtheta, -3*dtheta, -2*dtheta, -dtheta, 0.0, dtheta, 2*dtheta, 3*dtheta, 4*dtheta, 5*dtheta];
%contour_values = [1.0001 1.0011 1.0022 1.0022 1.0033 1.0044 1.0055 1.0065];
%contour_values = linspace(-0.01,0.01,41) / 288.15;
title_true = 1;

%kmin = 50;
%kmax = 53;
kmin = 0;
kmax = 601;
dk   = 1;

if strcmp(test_case, 'Equatorial-Long-Wave')
    scalefactor = 2.0;
    ncx = 300; 
    ncy = 20;  
    %ncx = 600; 
    %ncy = 80;  
    L   = 8*3000.0 * scalefactor;  % 
    x0  = 0.0*L;
    H   = 10.0;  %
    aspect = [L/H/3 1 1];
    velosc = 100;  % velocity unit of RKLM code
elseif strcmp(test_case, 'Internal-Wave-Long-Wave')
    scalefactor = 2.0;
    ncx = 300; 
    ncy = 20;  
    %ncx = 600; 
    %ncy = 80;  
    L   = 3000.0 * scalefactor;  % 
    x0  = 0.0*L;
    H   = 10.0;  %
    aspect = [L/H/3 1 1];
    velosc = 100;  % velocity unit of RKLM code
elseif strcmp(test_case, 'Internal-Wave-Strong-Strat')
    scalefactor = 1.0;
    ncx = 1200; 
    ncy = 40;  
    L   = 3000.0 * scalefactor;  % 
    x0  = 0.5*L;
    H   = 10.0;  %
    aspect = [80 1 1];
    velosc = 100;  % velocity unit of RKLM code
elseif strcmp(test_case, 'Skamarock-Klemp-Internal-Wave')
    scalefactor = 1.0;
    ncx = 300; 
    ncy = 10;
    ncz = 5;
    L   = 300.0 * scalefactor;  % [km] 
    x0  = 0.5*L;
    H   = 10.0;    % [km] 
    B   = 20.0;    % [km]
    aspect = [16 1 1];
    velosc = 100;  % [m/s] velocity unit of RKLM code
elseif strcmp(test_case, 'Rising-Bubble')
    ncx = 160;  
    ncy =  80;  
    L   = 2.0;  
    x0  = 1.0;
    H   = 1.0; 
    aspect = [1 1 1];
    velosc = 100;  % velocity unit of RKLM code
elseif strcmp(test_case, 'Travelling-Vortex')
    ncx = 64;  % 512; 256;
    ncy = 64;  % 512; 256;
    ncz =  8;  % 512; 256;
    L   = 1.0;  
    x0  = 0.5;
    H   = 1.0;
    B   = 1.0/8.0;
    aspect = [1 1 1];
    velosc = 100;  % velocity unit of RKLM code
elseif strcmp(test_case, 'Straka')
    ncx = 512;  
    ncy = 64;  
    L  = 52.2;  % 
    x0 = 0.0*L;
    H  = 6.4;  %
    aspect = [1 1 1];
    velosc = 100;  % velocity unit of RKLM code
elseif strcmp(test_case, 'Smolarkiewicz-Margolin-Breaking-Wave')
    ncx = 240;  
    ncy = 120;  
    L   = 120.0;  % 
    x0  =  60.0;
    H   =  60.0;  %
    aspect = [1 1 1];
    velosc = 100;  % velocity unit of RKLM code
elseif strcmp(test_case, 'Advection')
    ncx = 31;  
    ncy = 31;  
    ncz = 31;  
    L  = 2.0;  % 
    x0 = 0.0*L;
    H  = 2.0;
    B  = 2.0;%
    aspect = [1 1 1];
    velosc = 100;  % velocity unit of RKLM code
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
if ncz == 1
    nnz = 1;
    dumsz = 0;
else
    nnz = ncz+1;
    dumsz = 2;
end

modelfigstr = strcat('  (',modelstr,')');

rhoY_diff = 0;
rhoZ_diff = 0;
transp    = 0;

if strcmp(slice, 'full3D')
    kkmin = 1;
    dkk   = 2;
    kkmax = ncz;
else
    nslice = floor(ncz/2);
end

folderstring = strcat('/Users/rupert/Documents/Computation/RKLM_Reference/low_Mach_gravity_',modelstr);

% cell-centered fields
%varstr = 'rho'; folderstr = 'rho'; titlestr = 'rho'; ndummy = 2; arraysize = [ncx ncy ncz];
%varstr = 'p'; folderstr = 'p'; titlestr = 'p'; ndummy = 2; arraysize = [ncx ncy ncz];
%varstr = 'S'; folderstr = 'S'; titlestr = 'S'; ndummy = 2; arraysize = [ncx ncy ncz];
%varstr = 'rhoY';  folderstr = 'rhoY'; titlestr = 'rhoY'; ndummy = 2; arraysize = [ncx ncy ncz]; rhoY_diff = 1;
%varstr = 'drhoY';  folderstr = 'drhoY'; titlestr = 'drhoY'; ndummy = 2; arraysize = [ncx ncy ncz];
%varstr = 'Y';  folderstr = 'Y'; titlestr = '\theta'; ndummy = 2; arraysize = [ncx ncy ncz];
%varstr = 'dY';  folderstr = 'dY'; titlestr = 'd\theta'; ndummy = 2; arraysize = [ncx ncy ncz];
%varstr = 'buoy';  folderstr = 'buoy'; titlestr = 'buoy'; ndummy = 2; arraysize = [ncx ncy ncz];
%varstr = 'rhoZp';  folderstr = 'rhoZp'; titlestr = 'rhoZp'; ndummy = 2; arraysize = [ncx ncy ncz];
%varstr = 'rhoZB';  folderstr = 'rhoZB'; titlestr = 'rhoZB'; ndummy = 2; arraysize = [ncx ncy ncz];
%varstr = 'rhoZ';  folderstr = 'rhoZ'; titlestr = 'rhoZ'; ndummy = 2; arraysize = [ncx ncy ncz];
%varstr = 'Z';  folderstr = 'Z'; titlestr = 'Z'; ndummy = 2; arraysize = [ncx ncy ncz]; rhoZ_diff = 0;
varstr = 'u';  folderstr = 'u'; titlestr = 'u'; ndummy = 2; arraysize = [ncx ncy ncz]; symmetry = -1*symmetry;
%varstr = 'v';  folderstr = 'v'; titlestr = 'v'; ndummy = 2; arraysize = [ncx ncy ncz];
%varstr = 'w';  folderstr = 'w'; titlestr = 'w'; ndummy = 2; arraysize = [ncx ncy ncz];
%varstr = 'qv';  folderstr = 'qv'; titlestr = 'qv'; ndummy = 2; arraysize = [ncx ncy ncz];
%varstr = 'qc';  folderstr = 'qc'; titlestr = 'qc'; ndummy = 2; arraysize = [ncx ncy ncz];
%varstr = 'qr';  folderstr = 'qr'; titlestr = 'qr'; ndummy = 2; arraysize = [ncx ncy ncz];
%varstr = 'flux_rhou';  folderstr = 'fluxes'; titlestr = 'flux_rhou'; ndummy = 2; arraysize = [ncx+1 ncy ncz]; symmetry = -1*symmetry;
%varstr = 'flux_rhov';  folderstr = 'fluxes'; titlestr = 'flux_rhov'; ndummy = 2; arraysize = [ncy+1 ncz ncx]; symmetry = -1*symmetry;
%varstr = 'flux_rhov';  folderstr = 'fluxes'; titlestr = 'flux_rhov'; ndummy = 2; arraysize = [ncz+1 ncx ncy]; symmetry = -1*symmetry;

%varstr = 'p2_c';  folderstr = 'p2_c'; titlestr = '\pi'; ndummy = 2; arraysize = [ncx ncy ncz];
%varstr = 'dp2_c';  folderstr = 'dp2_c'; titlestr = 'd\pi'; ndummy = 2; arraysize = [ncx ncy ncz];
%varstr = 'dpdim';  folderstr = 'dpdime'; titlestr = 'dp [Pa]'; ndummy = 2; arraysize = [ncx ncy ncz];
%varstr = 'rhs_cells';  folderstr = 'rhs_cells'; titlestr = 'rhs_c'; ndummy = 2; arraysize = [ncx ncy ncz];

%varstr = 'p2_n';  folderstr = 'p2_nodes'; titlestr = '\pi_n';    ndummy = 2; arraysize = [nnx nny nnz];
%varstr = 'dp2_n';  folderstr = 'dp2_nodes'; titlestr = 'd\pi_n';    ndummy = 2; arraysize = [nnx nny nnz];
%varstr = 'rhs_nodes';  folderstr = 'rhs_nodes'; titlestr = 'rhs_n';    ndummy = 2; arraysize = [nnx nny nnz];

%varstr = 'advflux_x';  folderstr = 'advflux'; titlestr = 'advflux_x'; ndummy = 2; arraysize = [ncx+1 ncy ncz]; symmetry = -1*symmetry;
%varstr = 'advflux_y';  folderstr = 'advflux'; titlestr = 'advflux_y'; ndummy = 2; arraysize = [ncy+1 ncz ncx]; symmetry = -1*symmetry; transp = 1;
%varstr = 'advflux_z';  folderstr = 'advflux'; titlestr = 'advflux_z'; ndummy = 2; arraysize = [ncz+1 ncx ncy]; symmetry = -1*symmetry; transp = 1;



scrsz = get(0,'ScreenSize');
figure1 = figure('Position',[1 scrsz(4)/1.5 scrsz(3)/1.5 scrsz(4)/1.75]);
if title_true
    texthandle = annotation('textbox', [0.45 0.95 0.2 0.05]);
    set(texthandle, 'String', strcat(titlestr,modelfigstr,' 3D'), 'FontSize', 18, 'FontName', 'Helvetica', 'LineStyle', 'none');
end
%set(texthandle, 'String', titlestr, 'FontSize', 14, 'FontName', 'Optima');

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
    v = hdfread(filestr, '/Data-Set-2', 'Index', {[1  1  1],[1  1  1],[arraysize(1)+2*ndummy  arraysize(2)+2*ndummy  arraysize(3)+2*ndummy]});
    
    if strcmp(slice, 'full3D')

        [nx, ny, nz] = size(v);
        nx = nx - 4;
        ny = ny - 4;
        nz = nz - 4;
        
        dx = L/nx;
        dy = H/ny;
        dz = B/nz;
        
        x = linspace(0.5*dx-(nx/2)*dx,-0.5*dx+(nx/2)*dx,nx);
        y = linspace(0.5*dy,-0.5*dy+ny*dy,ny);
        z = linspace(0.5*dz-(nz/2)*dz,-0.5*dz+(nz/2)*dz,nz);
                
        th = v(3:1:nx+2, 3:1:ny+2, 3:1:nz+2);
                
        % Create contours
        %[X,Y,Z] = meshgrid(x,y,z);
        %contourslice(x,y,z,th,0.0,0.5,0.0);
        %contourslice(X,Y,Z,th,[],[0.33, 0.5, 0.66],[], [1.001, 1.002, 1.003, 1.004, 1.005, 1.006]);

        %isovalue = 1.001;
        %fvc = isocaps(th,isovalue);
        %fvc = isocaps(x,y,z,th,isovalue);
        %patch(fvc);

        % Create contours

        if filledcontours
            for kk = kkmin:dkk:kkmax
                contourf(x,y,transpose(th(:,:,kk)),15,'LineColor','auto');
                colormap Jet;
                colorbar('FontSize',14,'FontName','Helvetica');
                set(gca,'DataAspectRatio', aspect, 'FontSize',18,'FontName','Helvetica');
                pause;
            end
        else
            if fixed_contours
                contour(x,z,th,contour_values,'LineColor','k');
            else
                contour(x,z,th,10,'LineColor','k');
            end
            set(gca,'DataAspectRatio', aspect, 'FontSize',18,'FontName','Helvetica');
        end

        if title_true
            title(strcat(titlestr,kstr));
        end
        % Create xlabel
        xlabel('x','FontSize',18,'FontName','Helvetica');
        
        % Create ylabel
        ylabel('y','FontSize',18,'FontName','Helvetica');

        % Create ylabel
        zlabel('z','FontSize',18,'FontName','Helvetica');

        pause
    else
        
        % for now let's take it slice by slice
        if strcmp(slice, 'xy')
            velo = v(:,:,nslice);
        elseif strcmp(slice, 'zy')
            velo = transpose(reshape(v(nslice,:,:), [arraysize(2)+2*ndummy, arraysize(3)+2*ndummy]));
        end
        
        if strcmp(varstr, 'rhs_nodes')
            velo(1:3,:) = 0.0*velo(1:3,:);
            velo(:,1:3) = 0.0*velo(:,1:3);
            velo(end-2:end,:) = 0.0*velo(end-2:end,:);
            velo(:,end-2:end) = 0.0*velo(:,end-2:end);
        end
        
        % eigenmode test case
        [nx, nz] = size(velo);
        nx = nx - 4;
        nz = nz - 4;
        
        dx = L/nx;
        dz = H/nz;
        
        dth = 8.3333e-04;
        x = linspace(0.5*dx-(nx/2)*dx,-0.5*dx+(nx/2)*dx,nx);
        z = linspace(0.5*dz,-0.5*dz+nz*dz,nz);
        Yt = transpose(velo);
        th = Yt(3:1:nz+2, 3:1:nx+2);
        
        %figure(1);
        %contour(x,z,th, [1+dth 1+2*dth 1+3*dth 1+4*dth 1+5*dth 1+6*dth 1+7*dth 1+8*dth], 'k');
        %set(gca,'DataAspectRatio',[1 1 1]);
        
        
        % Create contour
        if filledcontours
            contourf(x,z,th,15,'LineColor','auto');
            %contourf(x,z,th,[1.01 1.05 1.1 1.15 1.2 1.24],'LineColor','auto');
            colormap Jet
            colorbar('FontSize',14,'FontName','Helvetica')
        else
            if fixed_contours
                contour(x,z,th,contour_values,'LineColor','k');
            else
                contour(x,z,th,10,'LineColor','k');
            end
        end
        set(gca,'DataAspectRatio', aspect, 'FontSize',18,'FontName','Helvetica');
        if title_true
            title(strcat(titlestr,kstr));
        end
        % Create xlabel
        xlabel('x','FontSize',18,'FontName','Helvetica');
        
        % Create ylabel
        ylabel('z','FontSize',18,'FontName','Helvetica');
        
        pause
    end
end


