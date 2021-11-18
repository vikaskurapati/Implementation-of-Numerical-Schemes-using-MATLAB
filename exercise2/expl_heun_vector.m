function y = expl_heun_vector(y_0, dt, t_end, f)
n = floor(t_end/dt) + 1;
m = length(y_0);
y = zeros(n,m);
y(1,:) = y_0;
i = 1;
t = 0;
while (t < t_end)
    gradient1 = f(t, y);
    mid_step = y(i,:) + dt*gradient1;
    gradient2 = f(t+dt, mid_step);
    gradient = 0.5*(gradient1 + gradient2);
    y(i+1,:) = y(i,:) + dt*gradient;
    t = t+dt;
    i = i+1;
end