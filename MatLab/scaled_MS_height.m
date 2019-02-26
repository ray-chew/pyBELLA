function [hs,hsp] = scaled_MS_height(q,xi,its)

F = @(z) q.*cos(xi+z);

y = 0.0*xi;

for k=1:1:its
    y = F(y);
end

hs  = y;
hsp = - q*sin(xi+y)./(1+q*sin(xi+y));
end

