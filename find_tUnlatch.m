function [timeUnlatch,Error,t,w] = find_tUnlatch(R,mus,muk,mLatch,mLoad,latchDescription,springDescription,latchStartConditions,fricFormInfo)
%FIND_TUNLATCH Summary of this function goes here
%   Detailed explanation goes here
Etol = eps;
E = 1;


xLatch = latchStartConditions(1);
vLatch = latchStartConditions(2);
ys = yLoad(R,xLatch);


w1 = latchDescription{1};
rest = latchDescription{2};
if strcmpi(w1,'constant_force') || strcmpi(w1,'constant_force ') || strcmpi(w1,'linear_motor')
   timeSpan = [0, 10.*sqrt(2.*R./(Fl(latchDescription,mus,R,xLatch,ys,latchStartConditions(2),0,0,0)./mLatch))];
elseif strcmpi(w1,'constant_velocity') || strcmpi(w1,'constant_velocity ')
   timeSpan = [0, R./str2num(rest)];
end



options = odeset('Events',@(t,w) unlatching_end(t,w,R,mus,muk,mLatch,mLoad,latchDescription,springDescription,fricFormInfo));
try    
    [t,w,te,we,ie] = ode45(@(t,w) latchODE(t,w,R,mus,muk,mLatch,latchDescription,springDescription,Inf,mLoad,fricFormInfo), [0,Inf], latchStartConditions, options);
    mask = imag(t) ~= 0;
    t(mask) = [];
    w(mask,:) = [];
    timeUnlatch = t(end);
    Error = eps;
    te = t(end);
catch ME
    timeUnlatch = 0;
    Error = 0;
    t = [];
    w = [];
end
end

function [value,isterminal,direction] = unlatching_end(t,w,R,mus,muk,mLatch,mLoad,latchDescription,springDescription,fricFormInfo)
y = yLoad(R,w(1));
temp = latchODE(t,w,R,mus,muk,mLatch,latchDescription,springDescription,Inf,mLoad,fricFormInfo);
d2x = temp(2);
d2y = d2yLoad(R,w(1),w(2),d2x);
Fspring = Fs(springDescription,y);
netF = d2y.*mLoad;

if Fspring < 0
   value = 0;
else
    value = netF-Fspring;
end
isterminal = 1;
direction = 0;
end

function [diff] = findRoot(N,a,b,c,fricInfo)
diff = ((c+a.*N)./b - Ff(N, fricInfo)).^2;

end