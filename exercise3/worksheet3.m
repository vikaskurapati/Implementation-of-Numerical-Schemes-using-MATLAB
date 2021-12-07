close all; clear all; clc;
figure_count=1;
i=1;

%% a)
figure(figure_count)
hold on
for dt = [1/2.0,1/4.0,1/8.0,1/16.0,1/32.0]
    t = 0:dt:5;
    y_eul = expl_euler(1, dt, 5, @gradient1a);
    error_expl(1,i) = err_cal(y_eul,exp(-7*t), dt, 5);
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
%% Printing errors to MATLAB console
disp("Explicit Euler Method");
VarNames = {'δt', '1/2', '1/4', '1/8', '1/16', '1/32'};
name = ["Error"; "Error Reduction"];
error_expl_red = ["-",error_expl(1,1)/error_expl(1,2),error_expl(1,2)/error_expl(1,3),error_expl(1,3)/error_expl(1,4),error_expl(1,4)/error_expl(1,5)];
Error = [error_expl(1,:);error_expl_red(1,:)];
Result = table(name,Error(:,1),Error(:,2),Error(:,3),Error(:,4),Error(:,5), 'VariableNames',VarNames);
disp(Result);

%% b) Testing newton
%y = newton(1,@Gtest,@dGtest);

%% c) Implementing implicit Euler Method
figure(figure_count)
i = 1;
hold on
for dt = [1/2.0,1/4.0,1/8.0,1/16.0,1/32.0]
    t = 0.0:dt:5.0;
    y_impl = impl_euler(1.0, dt, 5.0, @f,@df);
    error_impl(1,i) = err_cal(y_impl,exp(-7*t), dt, 5);
    i = i+1;
    plot(t, y_impl, 'DisplayName',strcat('dt = ', sprintf('%.3f', dt)));
end
%ylim([-1,1]);
title("Solution of Implicit Euler");
xlabel("t");
ylabel("x");
legend('show');
hold off
figure_count = figure_count + 1;
%% Printing errors to MATLAB console
disp("Implicit Euler Method");
error_impl_red = ["-",error_impl(1,1)/error_impl(1,2),error_impl(1,2)/error_impl(1,3),error_impl(1,3)/error_impl(1,4),error_impl(1,4)/error_impl(1,5)];
Error = [error_impl(1,:);error_impl_red(1,:)];
Result = table(name,Error(:,1),Error(:,2),Error(:,3),Error(:,4),Error(:,5), 'VariableNames',VarNames);
disp(Result);

%% Functions
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