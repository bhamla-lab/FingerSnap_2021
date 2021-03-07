function [dy] = dyLoad(R,x,dx)
%DYLOAD Takes in Latch properties, position of latch, and velocity of latch
%and outputs the speed of the load.
%   Derivative of function in yLoad
    n=2;
    dy = R*(x^(n-1)*dx/R^n)*(1-(x/R)^n)^(1/n-1);

end

