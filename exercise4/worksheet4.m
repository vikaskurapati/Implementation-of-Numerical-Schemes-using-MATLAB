%%
figure_count = 1;
for nx = [3,7,15,31,63]
    ny = nx;
    %% Direct Solution with Full matrix
    [A, b] = matrices(nx,ny);
    solution_mat = A\b;
    T = zeros(nx+2,ny+2);
    for i = 2:nx+1
        for j = 2:ny+1
            T(i,j)=solution_mat((i-2)*nx+j-1);
        end
    end
    %% Direct Solution with Sparse matrix
    A_s = sparse(A);
    solution_mat_sparse = A_s\b;
    T_s = zeros(nx+2,ny+2);
    for i = 2:nx+1
        for j = 2:ny+1
            T_s(i,j)=solution_mat_sparse((i-2)*nx+j-1);
        end
    end
    %% Iterative Solution with Gauss Seidel
    x = linspace(0,1,nx+2);
    y = linspace(0,1,ny+2);
    [X,Y] = meshgrid(x,y);
%% Plotting matrix solved solution
    figure(figure_count)
    surf(X,Y,T);
    zlim([0,1.1]);
    colormap("jet");
    colorbar;
    figure_count = figure_count + 1;

    figure(figure_count)
    contourf(X,Y,T,'ShowText','on');
    colormap("jet");
    colorbar;
    figure_count = figure_count + 1;
%% Plotting solutions of sparse solvers
    figure(figure_count)
    surf(X,Y,T_s);
    zlim([0,1.1]);
    colormap("jet");
    colorbar;
    figure_count = figure_count + 1;

    figure(figure_count)
    contourf(X,Y,T_s,'ShowText','on');
    colormap("jet");
    colorbar;
    figure_count = figure_count + 1;

end