function [y] = yLoad(R,x)
%YLATCH Takes in Latch geometry and latch position and outputs load
%position.
%   Equation derived from formula for a circle
    n=2;
    y = R*(1-(1-(x/R).^n).^(1./n));

end

