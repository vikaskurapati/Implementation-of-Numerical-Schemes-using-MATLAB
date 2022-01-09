function [A, b] = matrices(nx,ny)
    n = (nx)*(ny);
    hx = 1/(nx+1);
    hy = 1/(ny+1);
    k = -2*((1/(hx*hx))+(1/(hy*hy)));
    A = eye(n)*k;
    b = zeros(n,1);
    counter=1;
    % anal = zeros(n,1);
    for i = 1:nx
        for j = 1:ny
            %%do the calculations
             %if(i~=1 && j~=1 && i~=nx && j~=ny)
             if(j~=1)
                 A(counter,counter-1) = 1/(hx*hx);
             end
             if(j~=ny)
                 A(counter,counter+1) = 1/(hx*hx);
             end
             if(i~=1)
                 A(counter,counter-ny) = 1/(hy*hy);
             end
             if(i~=nx)
                 A(counter,counter+ny) = 1/(hy*hy);
             end
                 b(counter) = -2*pi*pi*sin(pi*(i)/(nx+1))*sin(pi*(j)/(ny+1));
                 % anal(counter) = sin(pi*(i)/(nx+1))*sin(pi*(j)/(ny+1));
            counter = counter + 1;
        end
    end
    %A = sparse(A);
    %b = sparse(b);
    %x = A\b;
    %norm(b-A*x)/(n^0.5)
end