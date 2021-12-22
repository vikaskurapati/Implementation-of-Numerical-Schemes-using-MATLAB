%% Creating empty vectors to store data
figure_count = 1;
runtime_direct_full = [];
storage_direct_full = [];
runtime_direct_sparse = [];
storage_direct_sparse = [];
runtime_iterative_gauss = [];
storage_iterative_gauss = [];
error_gauss = [];
for nx = [3,7,15,31,63]
    ny = nx;
    %% Direct Solution with Full matrix
    [A, b] = matrices(nx,ny);
    tic
    solution_mat_direct = A\b;
    runtime_direct_full = [runtime_direct_full, toc];
    storage_direct_full = [storage_direct_full, numel(A)+numel(b)+numel(solution_mat_direct)];
    T = zeros(nx+2,ny+2);
    for i = 2:nx+1
        for j = 2:ny+1
            T(i,j)=solution_mat_direct((i-2)*nx+j-1);
        end
    end
    %% Direct Solution with Sparse matrix
    A_s = sparse(A);
    tic
    solution_mat_sparse = A_s\b;
    runtime_direct_sparse = [runtime_direct_sparse, toc];
    storage_direct_sparse = [storage_direct_sparse, nnz(A_s) + numel(b)+numel(solution_mat_sparse)];
    T_s = zeros(nx+2,ny+2);
    for i = 2:nx+1
        for j = 2:ny+1
            T_s(i,j)=solution_mat_sparse((i-2)*nx+j-1);
        end
    end
    %% Iterative Solution with Gauss Seidel
    tic
    solution_mat_g_s = gauss_seidel(b,nx,ny);    
    runtime_iterative_gauss = [runtime_iterative_gauss,toc];
    storage_iterative_gauss = [storage_iterative_gauss, numel(solution_mat_g_s)+numel(b)];
    T_g_s = zeros(nx+2,ny+2);
    for i = 2:nx+1
        for j = 2:ny+1
            T_g_s(i,j)=solution_mat_g_s(i-1,j-1);
        end
    end
    x = linspace(0,1,nx+2);
    y = linspace(0,1,ny+2);
    [X,Y] = meshgrid(x,y);
    T_anal = sin(pi*X).*sin(pi*Y);
    error_gauss = [error_gauss, sqrt(sum(sum((T_anal-T_g_s).*(T_anal-T_g_s)))/(nx*ny))];
    
    %% Plotting Direct Solution with Full matrix
    figure(figure_count)
    surf(X,Y,T);
    zlim([0,1.1]);
    title(strcat("Surface plot for direct solution with full matrix with grid size ", num2str(nx)));
    colormap("jet");
    colorbar;
    figure_count = figure_count + 1;

    figure(figure_count)
    contourf(X,Y,T,'ShowText','on');
    title(strcat("Contour plot for direct solution with sparse matrix with grid size ", num2str(nx)));
    colormap("jet");
    colorbar;
    figure_count = figure_count + 1;

%% Plotting Direct Solution with Sparse matrix
    figure(figure_count)
    surf(X,Y,T_s);
    zlim([0,1.1]);
    title(strcat("Surface plot for direct solution with sparse matrix with grid size ", num2str(nx)));
    colormap("jet");
    colorbar;
    figure_count = figure_count + 1;

    figure(figure_count)
    contourf(X,Y,T_s,'ShowText','on');
    title(strcat("Contour plot for direct solution with sparse matrix with grid size ", num2str(nx)));
    colormap("jet");
    colorbar;
    figure_count = figure_count + 1;

%% Plotting Iterative Solution with Gauss Seidel
    figure(figure_count)
    surf(X,Y,T_g_s);
    zlim([0,1.1]);
    title(strcat("Surface plot for iterative solution with Gauss Seidel with grid size ", num2str(nx)));
    colormap("jet");
    colorbar;
    figure_count = figure_count + 1;

    figure(figure_count)
    contourf(X,Y,T_g_s,'ShowText','on');
    title(strcat("Contour plot for direct solution with sparse matrix with grid size ", num2str(nx)));
    colormap("jet");
    colorbar;
    figure_count = figure_count + 1;
end
%% Printing Runtime and Storage to MATLAB console using tables for Direct Solver with full matrix
disp("Direct Solver with full matrix:");
VarNames = {'Nx=Ny=', '7', '15', '31', '63'};
name = ["Runtime"; "Storage"];
direct_full = [runtime_direct_full(1,2:5);storage_direct_full(1,2:5)];
Result_full = table(name,direct_full(:,1),direct_full(:,2),direct_full(:,3),direct_full(:,4), 'VariableNames',VarNames);
disp(Result_full);

%% Printing Runtime and Storage to MATLAB console using tables for Direct Solver with sparse matrix
disp("Direct Solver with Sparse matrix:");
direct_sparse = [runtime_direct_sparse(1,2:5);storage_direct_sparse(1,2:5)];
Result_sparse = table(name,direct_sparse(:,1),direct_sparse(:,2),direct_sparse(:,3),direct_sparse(:,4), 'VariableNames',VarNames);
disp(Result_sparse);

%% Printing Runtime and Storage to MATLAB console using tables for Iterative Solver with Gauss-Seidel
disp("Iterative Solver with Gauss-Seidel:");
iterative_gauss = [runtime_iterative_gauss(1,2:5);storage_iterative_gauss(1,2:5)];
Result_gauss = table(name,iterative_gauss(:,1),iterative_gauss(:,2),iterative_gauss(:,3),iterative_gauss(:,4), 'VariableNames',VarNames);
disp(Result_gauss);

%% Adding Nx = Ny = 127 to gauss seidel
nx = 127; ny = 127;
[A, b] = matrices(nx,ny);
solution_mat_g_s = gauss_seidel(b,nx,ny);    
runtime_iterative_gauss = [runtime_iterative_gauss,toc];
storage_iterative_gauss = [runtime_iterative_gauss, numel(solution_mat_g_s)+numel(b)];
T_g_s = zeros(nx+2,ny+2);
for i = 2:nx+1
    for j = 2:ny+1
        T_g_s(i,j)=solution_mat_g_s(i-1,j-1);
    end
end
x = linspace(0,1,nx+2);
y = linspace(0,1,ny+2);
[X,Y] = meshgrid(x,y);
T_anal = sin(pi*X).*sin(pi*Y);
error_gauss = [error_gauss, sqrt(sum(sum((T_anal-T_g_s).*(T_anal-T_g_s)))/(nx*ny))];

%% Printing errors to MATLAB console using tables
disp("Gauss-Seidel Solver:");
VarNames = {'Nx=Ny=', '7', '15', '31', '63', '127'};
name = ["Error"; "Error Reduction"];
error_gauss_red = ["-",error_gauss(1,2)/error_gauss(1,3),error_gauss(1,3)/error_gauss(1,4),error_gauss(1,4)/error_gauss(1,5),error_gauss(1,5)/error_gauss(1,6)];
Error = [error_gauss(1,2:6);error_gauss_red(1,:)];
Result_error = table(name,Error(:,1),Error(:,2),Error(:,3),Error(:,4),Error(:,5), 'VariableNames',VarNames);
disp(Result_error);