% this m-file assumes Y to be the potential temperature
% distribution loaded directly from the .hdf-file for
% the Bryan-Bubble-output

slice = 'xy'; %  'xy' 'zy' 'full3D'
showmode = 1;
filledcontours = 1;
fixed_contours = 0;
%contour_values = [1.0001 1.0011 1.0022 1.0022 1.0033 1.0044 1.0055 1.0065];
contour_values = linspace(-0.01,0.01,41) / 288.15;
title_true = 1;

kmin = 1;
kmax = 11;
dk   = 1;
%modelstr = '';
%modelstr = 'psinc';
modelstr = 'comp';

if strcmp(slice, 'xy')
    ncx = 160;  nnx = ncx+1;
    ncy =  80;  nny = ncy+1;
    ncz =   1;  nnz = ncz+1;
elseif strcmp(slice, 'zy')
    ncx =   1;  nnx = ncx+1;
    ncy =  80;  nny = ncy+1;
    ncz = 160;  nnz = ncz+1;
elseif strcmp(slice, 'full3D')
    ncx =  80;  nnx = ncx+1;
    ncy =  40;  nny = ncy+1;
    ncz =  80;  nnz = ncz+1;
    
    kkmin = ncz/2;
    dkk   =  1;
    kkmax = ncz/2;
end

L = 2;  %
H = 1;  %
B = 2;  %

aspect = [1 1 1];
%aspect = [10 1 1];

modelfigstr = strcat('  (',modelstr,')');

% cell-centered fields
%varstr = 'rho'; folderstr = 'rho'; titlestr = 'rho'; ndummy = 2; arraysize = [ncx ncy ncz];
%varstr = 'rhoY';  folderstr = 'rhoY'; titlestr = 'rhoY'; ndummy = 2; arraysize = [ncx ncy ncz];
%varstr = 'rhoZ';  folderstr = 'rhoZ'; titlestr = 'rhoZ'; ndummy = 2; arraysize = [ncx ncy ncz];
varstr = 'Y';  folderstr = 'Y'; titlestr = '\theta'; ndummy = 2; arraysize = [ncx ncy ncz];
%varstr = 'dY';  folderstr = 'dY'; titlestr = 'd\theta'; ndummy = 2; arraysize = [ncx ncy ncz];
%varstr = 'Z';  folderstr = 'Z'; titlestr = 'Z'; ndummy = 2; arraysize = [ncx ncy ncz];
%varstr = 'u';  folderstr = 'u'; titlestr = 'u'; ndummy = 2; arraysize = [ncx ncy ncz];
%varstr = 'v';  folderstr = 'v'; titlestr = 'w'; ndummy = 2; arraysize = [ncx ncy ncz];
%varstr = 'qv';  folderstr = 'qv'; titlestr = 'qv'; ndummy = 2; arraysize = [ncx ncy ncz];
%varstr = 'qc';  folderstr = 'qc'; titlestr = 'qc'; ndummy = 2; arraysize = [ncx ncy ncz];
%varstr = 'qr';  folderstr = 'qr'; titlestr = 'qr'; ndummy = 2; arraysize = [ncx ncy ncz];
%varstr = 'dp2_c';  folderstr = 'dp2_c'; titlestr = 'd\pi'; ndummy = 2; arraysize = [ncx ncy ncz];
%varstr = 'p2_c';  folderstr = 'p2_c'; titlestr = 'p2_c'; ndummy = 2; arraysize = [ncx ncy ncz];
%varstr = 'dpdim';  folderstr = 'dpdime'; titlestr = 'dp [Pa]'; ndummy = 2; arraysize = [ncx ncy ncz];

% node-centered fields
%varstr = 'dp2_n';  folderstr = 'p2_nodes'; titlestr = 'd\pi';    ndummy = 2; arraysize = [nnx nny nnz];
%varstr = 'rhs_nodes';  folderstr = 'rhs_nodes'; titlestr = 'rhs_n';    ndummy = 2; arraysize = [nnx nny nnz];

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
        filestr = strcat('/Users/work/Documents/Computation/LowMach_MG/low_Mach_gravity_',modelstr,'/',folderstr,'/',varstr,'_00',kstr,'.hdf');
    else if k < 100
            filestr = strcat('/Users/work/Documents/Computation/LowMach_MG/low_Mach_gravity_',modelstr,'/',folderstr,'/',varstr,'_0',kstr,'.hdf');
        else
            filestr = strcat('/Users/work/Documents/Computation/LowMach_MG/low_Mach_gravity_',modelstr,'/',folderstr,'/',varstr,'_',kstr,'.hdf');
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
            end
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
        ylabel('y','FontSize',18,'FontName','Helvetica');

        % Create ylabel
        zlabel('z','FontSize',18,'FontName','Helvetica');

        pause
    else
        
        % for now let's take it slice by slice
        if strcmp(slice, 'xy')
            velo = v(:,:,ndummy+1);
        elseif strcmp(slice, 'zy')
            velo = transpose(reshape(v(ndummy+1,:,:), [arraysize(2)+2*ndummy, arraysize(3)+2*ndummy]));
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


