function T = gauss_seidel(b,nx,ny)
    %% 
    tol = 1e-4;
    error = 1e1;
    T = zeros(nx,ny);
    T_old = T;
    hx = 1/(nx+1);
    hy = 1/(ny+1);
    k = 2*((1/(hx*hx))+(1/(hy*hy)));
    while(error>tol)
        %% Gauss Seidal iteration
        for i = 1:nx
            for j = 1:ny
                if(i~=1 && j~=1 && i~=nx && j~=ny)
                   T(i,j) = (1/k)*((T(i-1,j)+T(i+1,j))/(hx*hx) + (T(i,j-1)+T(i,j+1))/(hy*hy) - b((i-1)*nx+j));
                end
                if(i==1)
                    if(j==1)
                        T(i,j) = (1/k)*(T(i+1,j)/(hx*hx) + T(i,j+1)/(hy*hy) - b((i-1)*nx+j));
                    elseif (j==ny)
                        T(i,j) = (1/k)*((T(i+1,j))/(hx*hx) + (T(i,j-1))/(hy*hy) - b((i-1)*nx+j));
                    else
                        T(i,j) = (1/k)*((T(i+1,j))/(hx*hx) + (T(i,j-1)+T(i,j+1))/(hy*hy) - b((i-1)*nx+j));
                    end
                end
                if(i==nx)
                    if(j==1)
                        T(i,j) = (1/k)*((T(i-1,j))/(hx*hx) + (T(i,j+1))/(hy*hy) - b((i-1)*nx+j));
                    elseif (j==ny)
                        T(i,j) = (1/k)*((T(i-1,j))/(hx*hx) + (T(i,j-1))/(hy*hy) - b((i-1)*nx+j));
                    else
                        T(i,j) = (1/k)*((T(i-1,j))/(hx*hx) + (T(i,j-1)+T(i,j+1))/(hy*hy) - b((i-1)*nx+j));
                    end
                end
                if(j==1 && i~=1 && i~=nx)
                        T(i,j) = (1/k)*((T(i-1,j)+T(i+1,j))/(hx*hx) + (T(i,j+1))/(hy*hy) - b((i-1)*nx+j));
                end
                if(j==ny && i~=1 && i~=nx)
                        T(i,j) = (1/k)*((T(i-1,j)+T(i+1,j))/(hx*hx) + (T(i,j-1))/(hy*hy) - b((i-1)*nx+j));
                end
            end
        end
        %% Error Calculation
        ax = zeros(nx*ny,1);
        for i = 1:nx
            for j = 1:ny
                if(i~=1 && j~=1 && i~=nx && j~=ny)
                   ax((i-1)*nx+j) = (T(i-1,j)+T(i+1,j))/(hx*hx) + (T(i,j-1)+T(i,j+1))/(hy*hy) - k*T(i,j);
                end
                if(i==1)
                    if(j==1)
                        ax((i-1)*nx+j) = T(i+1,j)/(hx*hx) + T(i,j+1)/(hy*hy) - k*T(i,j);
                    elseif (j==ny)
                        ax((i-1)*nx+j) = (T(i+1,j))/(hx*hx) + (T(i,j-1))/(hy*hy) - k*T(i,j);
                    else
                        ax((i-1)*nx+j) = (T(i+1,j))/(hx*hx) + (T(i,j-1)+T(i,j+1))/(hy*hy) - k*T(i,j);
                    end
                end
                if(i==nx)
                    if(j==1)
                        ax((i-1)*nx+j) = (T(i-1,j))/(hx*hx) + (T(i,j+1))/(hy*hy) - k*T(i,j);
                    elseif (j==ny)
                        ax((i-1)*nx+j) = (T(i-1,j))/(hx*hx) + (T(i,j-1))/(hy*hy) - k*T(i,j);
                    else
                        ax((i-1)*nx+j) = (T(i-1,j))/(hx*hx) + (T(i,j-1)+T(i,j+1))/(hy*hy) - k*T(i,j);
                    end
                end
                if(j==1 && i~=1 && i~=nx)
                        ax((i-1)*nx+j) = (T(i-1,j)+T(i+1,j))/(hx*hx) + (T(i,j+1))/(hy*hy) - k*T(i,j);
                end
                if(j==ny && i~=1 && i~=nx)
                        ax((i-1)*nx+j) = (T(i-1,j)+T(i+1,j))/(hx*hx) + (T(i,j-1))/(hy*hy) - k*T(i,j);
                end
            end
        end
        error = norm(b-ax)/sqrt(nx*ny);
    end
end