function SymmetryTests(v,n)

aspect = [1 1 1];

vt = transpose(v);

% symmetry in the first index
thsymy = 0.5 * (vt - vt(end:-1:1,:));
thasyy = 0.5 * (vt + vt(end:-1:1,:));

% symmetry in the second index
thsymx = 0.5 * (vt - vt(:,end:-1:1));
thasyx = 0.5 * (vt + vt(:,end:-1:1));

% mirror symmetry w.r.t. center of domain
thsymm = 0.5 * (vt - vt(end:-1:1,end:-1:1));
thasym = 0.5 * (vt + vt(end:-1:1,end:-1:1));


figure(n)
subplot(2,3,1)
contourf(thsymx,20,'LineColor','auto');
colorbar
colormap Jet
set(gca,'DataAspectRatio', aspect);
title('sym x')
%
subplot(2,3,4)
contourf(thasyx,20,'LineColor','auto');
colorbar
colormap Jet
set(gca,'DataAspectRatio', aspect);
title('asy x')
%
subplot(2,3,2)
contourf(thsymy,20,'LineColor','auto');
colorbar
colormap Jet
set(gca,'DataAspectRatio', aspect);
title('sym y')
%
subplot(2,3,5)
contourf(thasyy,20,'LineColor','auto');
colorbar
colormap Jet
set(gca,'DataAspectRatio', aspect);
title('asy y')
%
subplot(2,3,3)
contourf(thsymm,20,'LineColor','auto');
colorbar
colormap Jet
set(gca,'DataAspectRatio', aspect);
title('sym m')
%
subplot(2,3,6)
contourf(thasym,20,'LineColor','auto');
colorbar
colormap Jet
set(gca,'DataAspectRatio', aspect);
title('asy m')
