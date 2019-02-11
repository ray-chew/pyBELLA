n   = 1000;
its = 5;

% a0 = 1.0;
% om = 1;
% u0 = 1;
% L  = 1;
% kx = 2*pi/L;
% kz = 1;
% q  = sqrt(2*a0*om)/u0/sqrt(kx^2+ kz^2);

kz = 0.5;
kx = 1.0;
q  = 0.25;

x  = linspace(0,2*pi,n);
xi = kx*x;
[h,hp] = scaled_MS_height(q,xi,its);
h  = h/kz;
hp = kx*hp/kz;

figure(11)
subplot(2,1,1)
plot(x, h, 'LineWidth', 1);
axis tight;
ylabel('h(x)');
set(gca, 'FontSize',18,'FontName','Helvetica');
subplot(2,1,2)
plot(x, hp, 'LineWidth', 1);
axis tight;
xlabel('x');
ylabel('h^\prime(x)');
set(gca, 'FontSize',18,'FontName','Helvetica');