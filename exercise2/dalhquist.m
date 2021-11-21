function dalhquist()
t = 0:0.1:5;
figure(1);
plot(t, exp(-t), 'DisplayName','Analytical');
xlabel('t');
ylabel('-e^t');
title('Analytical solution');
    function y = gradient(t, x)
        y = -x;
    end
    function final_err = err_cal(y_cal,y_analytical, dt, tend)
        err = 0;
        n = length(y_cal);
        for j=1:n
            err = err + (y_cal(j)-y_analytical(j))*(y_cal(j)-y_analytical(j));
        end
        final_err = sqrt(dt*err/tend);
    end
fig_count = 2;
figure(fig_count)
error = zeros(3,4);
i = 1;
for dt = [1,1/2.0,1/4.0,1/8.0]
    t = 0:dt:5;
    y_eul = expl_euler(1, dt, 5, @gradient);
    error(1,i) = err_cal(y_eul,exp(t), dt, 5);
    i = i+1;
    hold on
    plot(t, y_eul, 'DisplayName',strcat('dt = ', sprintf('%.6f', dt)));
end
plot(t, exp(-t), 'DisplayName','Analytical');
legend('show');
fig_count = fig_count + 1;
hold off

figure(fig_count)
i=1;
for dt = [1,1/2.0,1/4.0,1/8.0]
    y_heun = expl_heun(1,dt,5, @gradient);
    t = 0:dt:5;
    error(2,i) = err_cal(y_heun,exp(t), dt, 5);
    i = i+1;
    hold on
    plot(t, y_heun, 'DisplayName',strcat('dt = ', sprintf('%.6f', dt)));
end
plot(t, exp(-t), 'DisplayName','Analytical');
legend('show');
fig_count = fig_count + 1;
hold off

figure(fig_count)
i = 1;
for dt = [1,1/2.0,1/4.0,1/8.0]
    y_runge_kutta = expl_runge_kutta(1, dt, 5, @gradient);
    t = 0:dt:5;
    error(3,i) = err_cal(y_runge_kutta,exp(t), dt, 5);
    i = i+1;
    hold on
    plot(t, y_runge_kutta, 'DisplayName',strcat('dt = ', sprintf('%.6f', dt)));
end
plot(t, exp(-t), 'DisplayName','Analytical');
legend('show');
fig_count = fig_count + 1;
hold off

error %error printing and reduction printing left

% figure(fig_count)
% dt = [1,1/2.0,1/4.0,1/8.0];
% hold on
% plot(log(dt),log(error(1,:)));
% plot(log(dt), log(error(2,:)));
% plot(log(dt), log(error(3,:)));
% legend('Euler', 'Heun', 'Runge-Kutta');
end