function [ t, y ] = DE2_lalania9(p, q, g, t0, tN, y0, y1, h)

N = floor((tN-t0)/h);

t = zeros(1, N+1);
y = zeros(1, length(t));
dy = zeros(1, length(t));
y(1) = y0;
y(2)  = y0 + y1*h;
t(1) = t0;
t(2) = t0 + h;

dy(2) = (y(2)-y(1))/h;

for n = 2:N
    % y(n+1) = 2*y(n) - y(n-1) + y''(n) * h^2 is determined using the
    % approximations provided
    y(n+1) = 2*y(n) - y(n-1) + (-p(t(n))*dy(n) - q(t(n))*y(n) + g(t(n)))*h^2; 
    t(n+1) = t(n) + h;
    dy(n+1) = (y(n+1) - y(n))/h; 
end

end