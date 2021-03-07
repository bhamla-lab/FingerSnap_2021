function [d2y] = d2yLoad(R,x,dx,d2x)
%D2YLOAD Takes in Latch properties, position of latch, velocity of the
%latch, and acceleration of the latch and outputs acceleration of load.
%   Derivative of equation in dyLoad.
    n=2;
    d2y = R*( (((n-1).*x^(n-2).*(dx)^2+x^(n-1)*d2x)*(1-(x/R)^n)^(1/n-1))/R^n+((x^(n-1)*dx)*(1/n-1)*(1-(x/R)^n)^(1/n-2)*(-n*x^(n-1)*dx))/((R^n)^2));

end

