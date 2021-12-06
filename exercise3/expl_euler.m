function y= expl_euler(y_0, dt, t_end, f)
n = floor(t_end/dt) + 1;
y = zeros(n,1);
y(1) = y_0;
i = 1;
t = 0;
while (t < t_end)
    gradient = f(t, y(i));
    y(i+1) = y(i) + dt*gradient;
    t = t+dt;
    i = i+1;
end