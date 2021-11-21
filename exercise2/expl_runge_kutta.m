function y = expl_runge_kutta(y_0, dt, t_end, f)
n = floor(t_end/dt) + 1;
y = zeros(n,1);
y(1) = y_0;
i = 1;
t = 0;
while (t < t_end)
    gradient1 = f(t, y(i));
    midstep1 = y(i) + 0.5*dt*gradient1;
    gradient2 = f(t+0.5*dt, midstep1);
    midstep2 = y(i) + 0.5*dt*gradient2;
    gradient3 = f(t+0.5*dt, midstep2);
    midstep3 = y(i) + dt*gradient3;
    gradient4 = f(t+0.5*dt, midstep3);
    gradient = (gradient1 + 2*gradient2 + 2*gradient3 + gradient4)/6.0;
    y(i+1) = y(i) + dt*gradient;
    t = t+dt;
    i = i+1;
end