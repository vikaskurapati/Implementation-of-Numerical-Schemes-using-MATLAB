close all; clear all; clc;
figure_count=1;
i=1;
%%1a
figure(figure_count)
hold on
for dt = [1/2.0,1/4.0,1/8.0,1/16.0,1/32.0]
    t = 0:dt:5;
    y_eul = expl_euler(1, dt, 5, @gradient1a);
    error(1,i) = err_cal(y_eul,exp(-7*t), dt, 5);
    i = i+1;
    plot(t, y_eul, 'DisplayName',strcat('dt = ', sprintf('%.3f', dt)));
end
legend('show');
ylim([-1,1]);
title("Solution of Explicit Euler");
xlabel("t");
ylabel("x");
hold off
figure_count = figure_count + 1;

%%Testing newton
%y = newton(1,@Gtest,@dGtest);
figure(figure_count)
i = 1;
hold on
for dt = [1/2.0,1/4.0,1/8.0,1/16.0,1/32.0]
    t = 0.0:dt:5.0;
    y_impl = impl_euler(1.0, dt, 5.0, @f,@df);
    %error(1,i) = err_cal(y_eul,exp(-7*t), dt, 5);
    i = i+1;
    plot(t, y_impl, 'DisplayName',strcat('dt = ', sprintf('%.3f', dt)));
end
ylim([-1,1]);
title("Solution of Implicit Euler");
xlabel("t");
ylabel("x");
legend('show');
hold off
figure_count = figure_count + 1;

function y = f(x)
    y = -7*x;
end

function y = df(x)
    y = -7;
end

function y = gradient1a(t, x)
    %defined this function to maintain consistency with the scheme
    %where the gradient can be a function of x and t
    y = -7*x;
end

function final_err = err_cal(y_cal,y_analytical, dt, tend)
        err = 0;
        n = length(y_cal);
        for j=1:n
            err = err + (y_cal(j)-y_analytical(j))*(y_cal(j)-y_analytical(j));
        end
        final_err = sqrt(dt*err/tend);
end

function y = Gtest(x)
    y = x*x*x - 3;
end
function y = dGtest(x)
    y = 3*x*x;
end