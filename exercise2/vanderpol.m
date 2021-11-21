function vanderpol()
    function grad = grad(t, y)
        grad = zeros(1,2);
        grad(1) = y(2);
        grad(2) = -y(1) + y(2)*(1-y(1)*y(1));
    end

    y_0 = [1,1];
    dt = 0.1;
    t_end = 20;
    t = 0:dt:t_end;
    y = expl_heun_vector(y_0, dt, t_end, @grad);
    hold all
    plot(t, y(:,1), '*', 'Color','r', 'DisplayName','y vs t');
    plot(t, y(:,2), '-', 'Color', 'b', 'DisplayName','x vs t');
    plot(y(:,1),y(:,2), 'o', 'Color','G', 'DisplayName', 'x vs y');
    legend('show');
end