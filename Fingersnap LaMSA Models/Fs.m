function [F] = Fs(springDescription,y)
%FS Returns spring force from description of spring and from position of
%load. 
%   First parses spring description to retrieve k (spring constant in N/m)
%   and yEQ (position of eqilibrium in m) and calculates spring force based 
%   on simple hookean dynamics. 
k = springDescription{1};
yEq = springDescription{2};
F = -k*(y-yEq);
end

