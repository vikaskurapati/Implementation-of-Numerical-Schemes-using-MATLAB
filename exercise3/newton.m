function y= newton(x_0,G,dG)
    i=0;
    tol=1;
    xcurr = x_0;
    while(i<100 && tol>1e-8)
      b =  dG(xcurr)\G(xcurr);
      xnext=xcurr-b;
      i=i+1;
      tol=norm((xcurr-xnext));
      xcurr=xnext;
    end
    y = xnext;
end

