y_0 = [0.01;0.01];

y = newton(y_0, @G, @dG);

y

function f = G(y)
    f = zeros(size(y));
    f(1) = y(2);
    f(2) = -y(1) + 4*(1-y(1)*y(1))*y(2);
end

function df = dG(y)
    df = zeros(length(y));
    df(1,:) = [0,1];
    df(2,:) = [-2*4*y(1)*y(2)-1, 4*(1-y(1)*y(1))];
end
