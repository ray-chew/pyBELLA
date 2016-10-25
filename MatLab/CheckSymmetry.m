function dmax = CheckSymmetry(th,x,z,no,symmetry)

delth = th;
lenvec = size(delth(1,:));
len   = lenvec(2);

for i=1:1:len
delth(:,i) = th(:,i)-symmetry*th(:,len+1-i);
end
scrsz = get(0,'ScreenSize');
contourf(x,z,delth,25,'LineColor','auto');
aspect = [8 1 1];
set(gca,'DataAspectRatio', aspect, 'FontSize',18,'FontName','Helvetica');
xlabel('x','FontSize',18,'FontName','Helvetica');
ylabel('z','FontSize',18,'FontName','Helvetica');
title(num2str(no));

dmax = max(max(delth));

end