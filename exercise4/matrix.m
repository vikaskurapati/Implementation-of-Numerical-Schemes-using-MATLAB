 function [A, b] = matrix(nx,ny)
    n = (nx+2)*(ny+2);
    hx = 1/(nx+1);
    hy = 1/(ny+1);
    k = -2*((1/(hx*hx))+(1/(hy*hy)));
    A = eye(n);
    b = zeros(n,1);
    counter=1;
    anal = zeros(n,1);
    for i = 1:nx+2
        for j = 1:ny+2
            %%do the calculations
             if(i~=1 && j~=1 && i~=nx+2 && j~=ny+2)
                 A(counter,counter) = k;
                 A(counter,counter-1) = 1/(hy*hy);
                 A(counter,counter+1) = 1/(hy*hy);
                 A(counter,counter-(nx+2)) = 1/(hx*hx);
                 A(counter,counter+(nx+2)) = 1/(hx*hx);
                 b(counter) = -2*pi*pi*sin(pi*(i-1)/(nx+1))*sin(pi*(j-1)/(ny+1));
                 anal(counter) = sin(pi*(i-1)/(nx+1))*sin(pi*(j-1)/(ny+1));
             end
            counter = counter + 1;
        end
    end
    %A = sparse(A);
    %b = sparse(b);
    %x = A\b;
    %norm(b-A*x)/(n^2)
    

end