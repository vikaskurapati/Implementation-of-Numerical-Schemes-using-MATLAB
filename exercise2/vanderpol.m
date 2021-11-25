function vanderpol()
    function grad = grad(t, y)
        %defined this function to maintain consistency with the scheme
        %where the gradient can be a function of x and t
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
    subplot(2,1,1);
    hold on
    plot(t, y(:,1), '*', 'Color','r', 'DisplayName','y vs t');
    plot(t, y(:,2), '-', 'Color', 'b', 'DisplayName','x vs t');
    hold off
    legend('show');
    subplot(2,1,2)
    hold on
    plot(y(:,1),y(:,2), 'o', 'Color','G', 'DisplayName', 'x vs y');
    hold off
    legend('show');
    sgtitle('Vander-Pol-Oscillator')
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
    title("Optional: Phase Space of y vs x")
    fig_count = fig_count + 1;
end