N = 1000;

rho0 = 0.5;

r = linspace(0.0,1,N);
dr = 1.0/(N-1);
p = 0.0*r;

gamma = 1.4;
Gamma = (gamma-1) / gamma;

dpdr =@(x) Gamma * (rho0 + 0.5*(1-x.^2).^6) .* (1-x).^12 .* x.^11;
%dpdr =@(x) Gamma * (rho0 + 0.0*(1-x.^2).^6) .* (1-x).^12 .* x.^11;
uth  =@(x) 1024*(1-x).^6 .* x.^6;

p(1) = 0.0;
for k = 2:N
    p(k) = p(k-1) + 0.5*dr* 1024*1024*(dpdr(r(k-1)) + dpdr(r(k)));
end

figure(50)
plot(r,p-p(N))

figure(52)
plot(r,uth(r))