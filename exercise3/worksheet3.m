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
title("Solution of Dahlquist's Equation: Explicit Euler");
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
title("Solution of Dahlquist's Equation: Implicit Euler (Newton)");
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

%% g) Vanderpol oscillator started
figure(figure_count)
y_0 = [1;1];
y_expl_vanderpol = expl_euler(y_0, 0.1,20,@gradientg);
y_expl_vanderpol_2 = expl_euler(y_0, 0.05,20, @gradientg);
dt = 0.05;
t_end = 20;
t = 0:dt:t_end;
hold on
plot(t',y_expl_vanderpol_2(:,1), 'DisplayName','x');
plot(t',y_expl_vanderpol_2(:,2), 'DisplayName','y');
title("Solution of Van-der-Pol-Oscillator Equation: Explicit Euler");
xlabel("t");
ylabel("x & y");
hold off
legend('show');
figure_count = figure_count+1;

%% i) Vanderpol oscillator - Implicit
figure(figure_count)
y_0 = [1;1];
y_impl_vanderpol = impl_euler(y_0, 0.1, 20,@f_vand, @df_vand);
dt = 0.1;
t_end = 20;
t = 0:dt:t_end;
plot(y_impl_vanderpol);
figure_count = figure_count+1;
subplot(2,1,1);
hold on
plot(t, y_impl_vanderpol(1,:), 'DisplayName','x');
plot(t, y_impl_vanderpol(2,:), 'DisplayName','y');
xlabel("t");
ylabel("x & y");
hold off
legend('show');
subplot(2,1,2)
hold on
plot(y_impl_vanderpol(1,:),y_impl_vanderpol(2,:));
xlabel("x");
ylabel("y");
hold off
sgtitle('Solution of Van-der-Pol-Oscillator Equation: Implicit Euler (Newton)')
figure_count = figure_count + 1;

%% Functions
function grad = gradientg(t, y)
    grad = zeros(size(y));
    grad(1) = y(2);
    grad(2) = -y(1) + 4.0*y(2)*(1-y(1)*y(1));
end

function y = f_vand(x)
    y(1,:) = x(2);
    y(2,:) = -x(1) + 4.0*x(2)*(1-x(1)*x(1));
end

function y = df_vand(x)
    y(1,1) = 0;
    y(1,2) = 1;
    y(2,1) = -2.0*4.0*x(1)*x(2)-1.0;
    y(2,2) = 4.0*(1-x(1)*x(1));
end

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