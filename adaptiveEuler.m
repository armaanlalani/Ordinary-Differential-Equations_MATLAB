function [t, y] = adaptiveEuler(f, t0, tf, y0, h)

tol = 1e-8;
t = [t0];
y = [y0];

while t(end) + h < tf
    B = y(end) + h * f(t(end), y(end));
    C = y(end) + h/2 * f(t(end), y(end));
    C = C + h/2 * f(t(end) + h/2, C);
    A = C - B;
    
    if abs(A) <= tol
        t = [t; t(end) + h];
        y = [y; B];
    end
    h = 0.9 * h * min(max(tol / abs(A), 0.3), 2);
end