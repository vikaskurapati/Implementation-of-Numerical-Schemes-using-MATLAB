function vanderpol()
    function grad = grad(t, y)
        grad = zeros(1,2);
        grad(1) = y(2);
        grad(2) = -y(1) + y(2)*(1-y(1)*y(1));
    end
    global fig_count;
    y_0 = [1,1];
    dt = 0.1;
    t_end = 20;
    t = 0:dt:t_end;
    y = expl_heun_vector(y_0, dt, t_end, @grad);
    figure(fig_count)
    hold on
    plot(t, y(:,1), '*', 'Color','r', 'DisplayName','y vs t');
    plot(t, y(:,2), '-', 'Color', 'b', 'DisplayName','x vs t');
    plot(y(:,1),y(:,2), 'o', 'Color','G', 'DisplayName', 'x vs y');
    title('Van-der-Pol-Oscillator(Method of Heun)');
    legend('show');
    hold off
    fig_count = fig_count + 1;
        
    dxdt = @(x, y) y;
    dydt = @(x, y) (1-x*x)*y - x;

    x = -1:0.1:1;
    y = -1:0.1:1;
    [X, Y] = meshgrid(x,y);
    dX = dxdt(X,Y);
    dY = dydt(X,Y);
    figure(fig_count)
    quiver(X, Y, dX,dY);
    axis tight

    fig_count = fig_count + 1;
end