function [F] = Fl(latchDescription,mu,R,x,y,dx,dy,dx2,dy2)
%FL Outputs latch force given a latch description.
%   Parses w input to determine method for calculating latch force, then
%   solves for it. Additional unlatching methods i.e. hill motor can be
%   coded here.
w1 = latchDescription{1};


if strcmpi(w1,'constant_force')
    rest = latchDescription{2};
    F = rest;
elseif strcmpi(w1, 'linear_motor')
    rest1 = latchDescription{2};
    max_force = rest1(1);
    velocity = rest1(2);
    range_of_motion = rest1(3);
    braking  = latchDescription{3};
    
    if braking
        F = max_force*(1-dx/velocity).*(abs(x)<=range_of_motion);
    else
        F = max((max_force*(1-dx/velocity)) .* (abs(x-rest1(4))<=range_of_motion), 0);
    end
end
end

