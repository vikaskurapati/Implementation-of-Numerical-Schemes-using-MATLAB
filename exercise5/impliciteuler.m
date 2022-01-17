figure_count =1;
%figure('units','normalized','outerposition',[0 0 1 1])
f1 = figure;
hold all;
sgtitle("Implicit Euler, t=1/8")
%figure('units','normalized','outerposition',[0 0 1 1])
f2 = figure;
sgtitle("Implicit Euler, t=2/8")
%figure('units','normalized','outerposition',[0 0 1 1])
f3 = figure;
sgtitle("Implicit Euler, t=3/8")
%figure('units','normalized','outerposition',[0 0 1 1])
f4 = figure;
sgtitle("Implicit Euler, t=4/8")
for n = [3,7,15,31]
    nx = n;ny=n;
    dt = 1/64;
    A = matrices(nx, ny);
    I = eye(size(A));
    T = ones(nx*ny, 1);
    x = linspace(0,1,nx+2);
    y = linspace(0,1,ny+2);
    [X,Y] = meshgrid(x,y);
    t = 0;
    solution = [];
    while(t<1/2)
        T = gauss_seidel(I - dt*A , T);
        t = t+dt;
        if(abs(t-1/8)<1e-8 ||abs(t-1/4)<1e-8||abs(t-3/8)<1e-8||abs(t-1/2)<1e-8)
            solution = [solution, T];
        end
    end
    figure(f1)
    subplot(2,2,figure_count)
    plot_surface(X,Y,solution(:,1))
    title(strcat("N=",num2str(n),"& dt=",num2str(dt)))

    figure(f2)
    subplot(2,2,figure_count)
    plot_surface(X,Y,solution(:,2))
    title(strcat("N=",num2str(n),"& dt=",num2str(dt)))

    figure(f3)
    subplot(2,2,figure_count)
    plot_surface(X,Y,solution(:,3))
    title(strcat("N=",num2str(n),"& dt=",num2str(dt)))

    figure(f4)
    subplot(2,2,figure_count)
    plot_surface(X,Y,solution(:,4))
    title(strcat("N=",num2str(n),"& dt=",num2str(dt)))

    figure_count = figure_count +1;

end