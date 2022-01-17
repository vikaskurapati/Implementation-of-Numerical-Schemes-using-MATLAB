function A = matrices(nx,ny)
    n = (nx)*(ny);
    hx = 1/(nx+1);
    hy = 1/(ny+1);
    k = -2*((1/(hx*hx))+(1/(hy*hy)));
    A = eye(n)*k;
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
            counter = counter + 1;
        end
    end
end