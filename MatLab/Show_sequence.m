% this m-file assumes Y to be the potential temperature
% distribution loaded directly from the .hdf-file for
% the Bryan-Bubble-output

test_case = 'Skamarock-Klemp-Internal-Wave';
%test_case = 'Rising-Bubble';
%test_case = 'Smolarkiewicz-Margolin-Breaking-Wave';
%test_case = 'Straka';
%test_case = 'Travelling-Vortex';

showmode = 1;
filledcontours = 1;
fixed_contours = 0;
no_of_contours = 25;
symmetry = 0;        % in {0,1}
symmetrytest = 0;
showdummycells = 1;

th0 = -0.0015/300;
dth = 5e-4/300;
contour_values = [th0 th0+dth th0+2*dth th0+3*dth th0+4*dth th0+5*dth th0+6*dth th0+7*dth th0+8*dth th0+9*dth th0+10*dth];
%contour_values = [1.0001 1.0011 1.0022 1.0022 1.0033 1.0044 1.0055 1.0065];
%contour_values = linspace(-0.01,0.01,41) / 288.15;
title_true = 1;

kmin = 0;
kmax = 40;
dk   = 1;

%modelstr = '';
modelstr = 'psinc';
%modelstr = 'comp';


if strcmp(test_case, 'Skamarock-Klemp-Internal-Wave')
    scalefactor = 1.0;
    ncx = 600; 
    ncy = 40;  
    L   = 300.0 * scalefactor;  % 
    x0  = 0.5*L;
    H   = 10.0;  %
    aspect = [8 1 1];
elseif strcmp(test_case, 'Rising-Bubble')
    ncx = 160;  
    ncy =  80;  
    L   = 2.0;  
    x0  = 1.0;
    H   = 1.0; 
    aspect = [1 1 1];
elseif strcmp(test_case, 'Travelling-Vortex')
    ncx = 128;  % 512; 256;
    ncy = 128;  % 512; 256;
    L   = 1.0;  
    x0  = 0.5;
    H   = 1.0; 
    aspect = [1 1 1];
elseif strcmp(test_case, 'Straka')
    ncx = 512;  
    ncy = 64;  
    L  = 52.2;  % 
    x0 = 0.0*L;
    H  = 6.4;  %
    aspect = [1 1 1];
elseif strcmp(test_case, 'Smolarkiewicz-Margolin-Breaking-Wave')
    ncx = 240;  
    ncy = 120;  
    L   = 120.0;  % 
    x0  =  60.0;
    H   =  60.0;  %
    aspect = [1 1 1];
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

% cell-centered fields
%varstr = 'rho'; folderstr = 'rho'; titlestr = 'rho'; ndummy = 2; arraysize = [ncx ncy];
%varstr = 'p'; folderstr = 'p'; titlestr = 'p'; ndummy = 2; arraysize = [ncx ncy];
%varstr = 'S'; folderstr = 'S'; titlestr = 'S'; ndummy = 2; arraysize = [ncx ncy];
%varstr = 'rhoY';  folderstr = 'rhoY'; titlestr = 'rhoY'; ndummy = 2; arraysize = [ncx ncy]; rhoY_diff = 1;
%varstr = 'drhoY';  folderstr = 'drhoY'; titlestr = 'drhoY'; ndummy = 2; arraysize = [ncx ncy];
%varstr = 'Y';  folderstr = 'Y'; titlestr = '\theta'; ndummy = 2; arraysize = [ncx ncy];
%varstr = 'dY';  folderstr = 'dY'; titlestr = 'd\theta'; ndummy = 2; arraysize = [ncx ncy];
%varstr = 'rhoZp';  folderstr = 'rhoZp'; titlestr = 'rhoZp'; ndummy = 2; arraysize = [ncx ncy];
%varstr = 'rhoZB';  folderstr = 'rhoZB'; titlestr = 'rhoZB'; ndummy = 2; arraysize = [ncx ncy];
%varstr = 'rhoZ';  folderstr = 'rhoZ'; titlestr = 'rhoZ'; ndummy = 2; arraysize = [ncx ncy];
%varstr = 'Z';  folderstr = 'Z'; titlestr = 'Z'; ndummy = 2; arraysize = [ncx ncy]; rhoZ_diff = 0;
%varstr = 'u';  folderstr = 'u'; titlestr = 'u'; ndummy = 2; arraysize = [ncx ncy]; symmetry = -1*symmetry;
varstr = 'v';  folderstr = 'v'; titlestr = 'v'; ndummy = 2; arraysize = [ncx ncy];
%varstr = 'w';  folderstr = 'w'; titlestr = 'w'; ndummy = 2; arraysize = [ncx ncy];
%varstr = 'qv';  folderstr = 'qv'; titlestr = 'qv'; ndummy = 2; arraysize = [ncx ncy];
%varstr = 'qc';  folderstr = 'qc'; titlestr = 'qc'; ndummy = 2; arraysize = [ncx ncy];
%varstr = 'qr';  folderstr = 'qr'; titlestr = 'qr'; ndummy = 2; arraysize = [ncx ncy];
%varstr = 'p2_c';  folderstr = 'p2_c'; titlestr = '\pi'; ndummy = 2; arraysize = [ncx ncy];
%varstr = 'dp2_c';  folderstr = 'dp2_c'; titlestr = 'd\pi'; ndummy = 2; arraysize = [ncx ncy];
%varstr = 'dpdim';  folderstr = 'dpdime'; titlestr = 'dp [Pa]'; ndummy = 2; arraysize = [ncx ncy];
%varstr = 'rhs_cells';  folderstr = 'rhs_cells'; titlestr = 'rhs_c'; ndummy = 2; arraysize = [ncx ncy];

%varstr = 'dp2_n';  folderstr = 'dp2_nodes'; titlestr = 'd\pi_n';    ndummy = 2; arraysize = [nnx nny];
%varstr = 'p2_n';  folderstr = 'p2_nodes'; titlestr = '\pi_n';    ndummy = 2; arraysize = [nnx nny];
%varstr = 'rhs_nodes';  folderstr = 'rhs_nodes'; titlestr = 'rhs_n';    ndummy = 2; arraysize = [nnx nny];

%varstr = 'advflux_x';  folderstr = 'advflux'; titlestr = 'advflux_x'; ndummy = 2; arraysize = [ncx+1 ncy]; symmetry = -1*symmetry;
%varstr = 'advflux_y';  folderstr = 'advflux'; titlestr = 'advflux_y'; ndummy = 2; arraysize = [ncy+1 ncx]; symmetry = -1*symmetry; transp = 1;


scrsz = get(0,'ScreenSize');
figure1 = figure('Position',[1 scrsz(4)/2 scrsz(3)/1 scrsz(4)/2.5]);
if abs(symmetry) == 1 
    figure2 = figure('Position',[scrsz(4)/2 0 scrsz(3)/1 scrsz(4)/2.5]);
end

if title_true == 0
    texthandle = annotation('textbox', [0.45 0.95 0.2 0.05]);
    set(texthandle, 'String', strcat(titlestr,modelfigstr,' 2D'), 'FontSize', 18, 'FontName', 'Helvetica', 'LineStyle', 'none');
end
%set(texthandle, 'String', titlestr, 'FontSize', 14, 'FontName', 'Optima');

for k = kmin:dk:kmax
    kstr = num2str(k);
    if k < 10
        filestr = strcat(folderstring,'/',folderstr,'/',varstr,'_00',kstr,'.hdf');
    else if k < 100
            filestr = strcat(folderstring,'/',folderstr,'/',varstr,'_0',kstr,'.hdf');
        else
            filestr = strcat(folderstring,'/',folderstr,'/',varstr,'_',kstr,'.hdf');
        end
    end
    if nny == 1
        v = hdfread(filestr, '/Data-Set-2', 'Index', {1,1,arraysize(1)+dumsx*ndummy});
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
        xlabel('x','FontSize',18,'FontName','Helvetica');
        
        % Create ylabel
        ylabel(titlestr,'FontSize',18,'FontName','Helvetica');

    else
        v = hdfread(filestr, '/Data-Set-2', 'Index', {[1  1],[1  1],[arraysize(1)+dumsx*ndummy  arraysize(2)+dumsy*ndummy]});
        
        if transp == 1
            v = transpose(v);
        end

        if rhoY_diff
            velo = v - rhoY_diff;
        elseif rhoZ_diff
            velo = v - min(min(v));
        else
            velo = v;
        end
%         if strcmp(varstr, 'rhs_nodes')
%             velo(1:3,:) = 0.0*velo(1:3,:);
%             velo(:,1:3) = 0.0*velo(:,1:3);
%             velo(end-2:end,:) = 0.0*velo(end-2:end,:);
%             velo(:,end-2:end) = 0.0*velo(:,end-2:end);
%         end
        
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
            contourf(x,z,th,25,'LineColor','auto');
            colormap Jet
            colorbar('FontSize',14,'FontName','Helvetica')
        else
            if fixed_contours
                contour(x,z,th,contour_values,'LineColor','k');
            else
                contour(x,z,th,no_of_contours,'LineColor','k');
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
        
        if abs(symmetry) == 1
            figure(figure2);
            NonSymmAmpl = CheckSymmetry(th,x,z,k,symmetry);
            colormap Jet
            colorbar('FontSize',14,'FontName','Helvetica')
            %NonSymmAmpl;
        end

    end
         
    if symmetrytest == 1
        SymmetryTests(transpose(th),55);
        figure(figure1)
    end
    pause
end


