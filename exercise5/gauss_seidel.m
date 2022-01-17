function x = gauss_seidel(A, b)
    x = zeros(size(b));
    error = 1;
    for iter=1:100
        for j = 1:size(A,1)
            x(j) = (b(j) - sum(A(j,:)'.*x) + A(j,j)*x(j))/A(j,j);
        end
        error = norm(A*x - b, 2);
        if error > 1e-6
            break;
        end
    end
end