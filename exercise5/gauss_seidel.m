function x = gauss_seidel(A, b)
    x = b;
    error = 1;
    for iter=1:100
        for i = 1:size(A,1)
            x(i) = (b(i) - dot(A(i,:),x) + A(i,i)*x(i))/A(i,i);
        end
        error = norm(A*x - b, 2);
        if error < 1e-6
            break;
        end
    end
end