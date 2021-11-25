function dalhquist()
t = 0:0.1:5;
global fig_count;
figure(fig_count);
plot(t, exp(-t), 'DisplayName','Analytical');
xlabel('t');
ylabel('e^{-t}');
title('Analytical solution');
    function y = gradient(t, x)
        %defined this function to maintain consistency with the scheme
        %where the gradient can be a function of x and t
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
fig_count = fig_count+1;
figure(fig_count)
error = zeros(3,4);
i = 1;
for dt = [1,1/2.0,1/4.0,1/8.0]
    t = 0:dt:5;
    y_eul = expl_euler(1, dt, 5, @gradient);
    error(1,i) = err_cal(y_eul,exp(-t), dt, 5);
    i = i+1;
    hold on
    plot(t, y_eul, 'DisplayName',strcat('dt = ', sprintf('%.3f', dt)));
end
plot(t, exp(-t), 'DisplayName','Analytical');
xlabel('t');
ylabel('x');
title('Dahlquists test Equation: Explicit Euler method');
legend('show');
fig_count = fig_count + 1;
hold off

figure(fig_count)
i=1;
for dt = [1,1/2.0,1/4.0,1/8.0]
    %Replaced the Heun with the extended heun for vectors to show that it
    %Works with scalars also
    y_heun = expl_heun_vector(1,dt,5, @gradient);
    t = 0:dt:5;
    error(2,i) = err_cal(y_heun,exp(-t), dt, 5);
    i = i+1;
    hold on
    plot(t, y_heun, 'DisplayName',strcat('dt = ', sprintf('%.3f', dt)));
end
plot(t, exp(-t), 'DisplayName','Analytical');
xlabel('t');
ylabel('x');
title('Dahlquists test Equation: Method of Heun');
legend('show');
fig_count = fig_count + 1;
hold off

figure(fig_count)
i = 1;
for dt = [1,1/2.0,1/4.0,1/8.0]
    y_runge_kutta = expl_runge_kutta(1, dt, 5, @gradient);
    t = 0:dt:5;
    error(3,i) = err_cal(y_runge_kutta,exp(-t), dt, 5);
    i = i+1;
    hold on
    plot(t, y_runge_kutta, 'DisplayName',strcat('dt = ', sprintf('%.3f', dt)));
end
plot(t, exp(-t), 'DisplayName','Analytical');
xlabel('t');
ylabel('x');
title('Dahlquists test Equation: Runge-Kutta method');
legend('show');
fig_count = fig_count + 1;
hold off

%% Error
disp('Explicit Euler method (q=1):');
disp('δt =      1         1/2         1/4          1/8 ');
fprintf('Error=  %.6f   %.6f    %.6f    %.6f\n',error(1,1),error(1,2),error(1,3),error(1,4));
fprintf('Error Red=   --    %.6f    %.6f    %.6f\n\n',error(1,2)/error(1,1),error(1,3)/error(1,2),error(1,4)/error(1,3));

disp('Method of Heun (q=2):');
disp('δt =      1         1/2         1/4          1/8 ');
fprintf('Error=  %.6f   %.6f    %.6f    %.6f\n',error(2,1),error(2,2),error(2,3),error(2,4));
fprintf('Error Red=   --    %.6f    %.6f    %.6f\n\n',error(2,2)/error(2,1),error(2,3)/error(2,2),error(2,4)/error(2,3));

disp('Runge-Kutta method (q=4):');
disp('δt =      1         1/2         1/4          1/8 ');
fprintf('Error=  %.6f   %.6f    %.6f    %.6f\n',error(3,1),error(3,2),error(3,3),error(3,4));
fprintf('Error Red=   --    %.6f    %.6f    %.6f\n',error(3,2)/error(3,1),error(3,3)/error(3,2),error(3,4)/error(3,3));

%% Log Log Plot for errors and dt to check the order of the schemes
% Uncomment the below lines to get the plots
% figure(fig_count)
% fig_count = fig_count+1;
% dt = [1,1/2.0,1/4.0,1/8.0];
% hold on
% plot(log(dt),log(error(1,:)));
% plot(log(dt), log(error(2,:)));
% plot(log(dt), log(error(3,:)));
% xlabel('log(dt)');
% ylabel('log(error)');
% title('Error Plot in Log Scale');
% legend('Euler', 'Heun', 'Runge-Kutta');
end