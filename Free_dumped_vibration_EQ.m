t = linspace(0,20,1000);
k = 20;
beta = 1;
m = 2.5;

w_n = sqrt(k/m);
zeta = beta/(2*w_n*m);

tol = 1e-14;
x_lim = 1e-4;   % stop when |x| becomes smaller than this value

if zeta < 1 - tol
    w_d = w_n*sqrt(1 - zeta^2);

    x = exp(-zeta*w_n.*t) .* sin(w_d.*t) ./ w_d;
    dx_x = exp(-zeta*w_n.*t) .* ...
           (cos(w_d.*t) - (zeta*w_n./w_d).*sin(w_d.*t));

elseif abs(zeta - 1) <= tol
    x = t .* exp(-w_n.*t);
    dx_x = exp(-w_n.*t) .* (1 - w_n.*t);

else
    s = sqrt(zeta^2 - 1);
    r1 = -w_n*(zeta - s);
    r2 = -w_n*(zeta + s);

    x = (exp(r1.*t) - exp(r2.*t)) ./ (r1 - r2);
    dx_x = (r1*exp(r1.*t) - r2*exp(r2.*t)) ./ (r1 - r2);
end

% Find local maxima
imax = find(dx_x(1:end-1) > 0 & dx_x(2:end) < 0);
a = length(imax);

% Stop only after the first maximum, when |x| becomes small enough
if ~isempty(imax)
    idx_after_peak = find(abs(x(imax(1):end)) <= x_lim, 1, 'first');

    if ~isempty(idx_after_peak)
        idx_stop = imax(1) + idx_after_peak - 1;
        t = t(1:idx_stop);
        x = x(1:idx_stop);
        dx_x = dx_x(1:idx_stop);

        imax = find(dx_x(1:end-1) > 0 & dx_x(2:end) < 0);
        a = length(imax);
    end
end

figure;
plot(t, x, 'LineWidth', 1.5);
hold on;
plot(t(imax), x(imax), 'ro', 'MarkerSize', 7, 'LineWidth', 1.5);
grid on;
xlabel('Time t [s]');
ylabel('Position x [m]');
title('Free damped vibration graph');