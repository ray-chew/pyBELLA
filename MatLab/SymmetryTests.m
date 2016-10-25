function SymmetryTests(v,n)

vt = transpose(v);

thsym = 0.5 * (vt - vt(:,end:-1:1));
thasy = 0.5 * (vt + vt(:,end:-1:1));

figure(n)
subplot(2,1,1)
contourf(thsym,20,'LineColor','auto');
colorbar
subplot(2,1,2)
contourf(thasy,20,'LineColor','auto');
colorbar