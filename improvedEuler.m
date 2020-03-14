function [t, y] = improvedEuler(f, t0, tf, y0, h)

t = t0 : h : tf;
N = length(t);
y = NaN(N, 1);

y(1) = y0;
for i = 1 : N - 1
    approx = y(i) + h * f(t(i), y(i));
    y(i + 1) = y(i) + h / 2 * (f(t(i), y(i)) + f(t(i + 1), approx));
end