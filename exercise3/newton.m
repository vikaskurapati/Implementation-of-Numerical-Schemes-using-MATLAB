function y= newton(x_0,G,dG)
    i=0;
    tol=1;
    xcurr = x_0;
    n = length(x_0);
    while(i<100 && tol>1e-8)
      if (rank(dG(xcurr)) == n)
        xnext=xcurr-dG(xcurr)\G(xcurr);
        i=i+1;
        tol=norm((xcurr-xnext)); %Using norm of the error vector as a tolerance limit
        xcurr=xnext;
      else
          xnext = ones(n,1)*NaN;
          break
      end
    end
    if (i == 100)
        disp('Newtons Method didnt converge, taking the last value from the iteration');
    end
    y = xnext;
end