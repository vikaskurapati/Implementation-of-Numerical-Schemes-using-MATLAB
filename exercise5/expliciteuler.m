figure_count =1;
figure('units','normalized','outerposition',[0 0 1 1])
hold all;
figure(1)
sgtitle("Explicit Euler,t=1/8")
figure('units','normalized','outerposition',[0 0 1 1])
figure(2)
sgtitle("Explicit Euler,t=2/8")
figure('units','normalized','outerposition',[0 0 1 1])
figure(3)
sgtitle("Explicit Euler,t=3/8")
figure('units','normalized','outerposition',[0 0 1 1])
figure(4)
sgtitle("Explicit Euler,t=4/8")
for n = [3,7,15,31]
    nx = n;
    ny = n;
    A = sparse(matrices(nx, ny));
%for n = [31]
    for dt = [1/64, 1/128, 1/256, 1/512, 1/1024, 1/2048, 1/4096]
    %for dt = [1/1024, 1/2048]
        T = ones(nx*ny, 1);
        t = 0;
        x = linspace(0,1,nx+2);
        y = linspace(0,1,ny+2);
        [X,Y] = meshgrid(x,y);
        solution = [];
        while (t < 1/2)
            T = T + dt*A*T;
            t = t + dt;
            if(abs(t-1/8)<1e-8 ||abs(t-1/4)<1e-8||abs(t-3/8)<1e-8||abs(t-1/2)<1e-8)
                solution = [solution, T];
                h1 = figure('Visible','off');
                plot_surface(X,Y,T)
                title(strcat("N=",num2str(n),"& dt=",num2str(dt),"& t=",num2str(t)))
                saveas(gca,strcat("N=",num2str(n),"& dt=",num2str(dt),"& t=",num2str(t),".png"))
            end
        end
        figure(1)
        subplot(4,7,figure_count)
        plot_surface(X,Y,solution(:,1))
        title(strcat("N=",num2str(n),"& dt=",num2str(dt)))

        figure(2)
        subplot(4,7,figure_count)
        plot_surface(X,Y,solution(:,2))
        title(strcat("N=",num2str(n),"& dt=",num2str(dt)))

        figure(3)
        subplot(4,7,figure_count)
        plot_surface(X,Y,solution(:,3))
        title(strcat("N=",num2str(n),"& dt=",num2str(dt)))

        figure(4)
        subplot(4,7,figure_count)
        plot_surface(X,Y,solution(:,4))
        title(strcat("N=",num2str(n),"& dt=",num2str(dt)))

        figure_count = figure_count +1;
    end
end