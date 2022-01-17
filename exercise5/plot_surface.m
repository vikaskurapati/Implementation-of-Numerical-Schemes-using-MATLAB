function plot_surface(X,Y,T)
T_s = zeros(size(X));

for i = 2:size(X,1)-1
    for j = 2:size(X, 2)-1
        T_s(i,j)=T((i-2)*(size(X,1)-2)+j-1);
    end
end
surf(X,Y,T_s)
end