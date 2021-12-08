function y = impl_euler(y_0, dt, t_end, f, df)
    t = 0.0:dt:t_end;
    %y = zeros(length(t),1);
    y(:,1) = y_0;
    I = eye(size(y_0,1));
    for i = 2:length(t)
        G = @(x) x - y(:,i-1) - dt*f(x);
        dG = @(x) I - dt*df(x);
        y(:,i) = newton(y(:,i-1), G, dG);
    end
end