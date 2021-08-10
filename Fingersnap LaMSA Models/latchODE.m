function [dw] = latchODE(t,w,R,mus,muk,mLatch,latchDescription,springDescription,timeUnlatch,mLoad,fricFormInfo)
%LATCHODE w = 2x1 matrix containing latch position and latch velocity.
%Takes w, R, and a description of latch kinematics (See Fl for
%more info) and outputs dw, a 2x1 matrix containing latch velocity and
%latch acceleration.
%   Based on differential equation derived from Newton's second law on
%   latch and load, and from geometry of the latch. Takes into account 
%   latch kinematics and force of friction. Calls functions Fl (describes 
%   latch kinematics) and Ff (horizontal component of friction on the latch)


tol = 1E-10;
options = optimset('Display','off');
y = yLoad(R,w(1));
[w1, rest] = strtok(latchDescription);

if strcmpi(w1,'constant_velocity')
    dw = [w(2);0];
else
    xLatch = w(1);
    vLatch = w(2);
    vLoad = dyLoad(R,xLatch,vLatch);
    sinTheta = xLatch./(sqrt(xLatch^2+(R-y).^2));
    cosTheta = (R-y)./(sqrt(xLatch^2+(R-y).^2)); 
    tanTheta = sinTheta./cosTheta;
    fricInfo = [mus, muk, xLatch, vLatch, y, vLoad, R, fricFormInfo];
    aLatch = 0;
    if t <= timeUnlatch && isreal(y)
       a = cosTheta + mLoad.*xLatch.*sinTheta./(mLatch.*R.*sqrt(1-(xLatch./R).^2));
       b = sinTheta + mLoad.*xLatch.*cosTheta./(mLatch.*R.*sqrt(1-(xLatch./R).^2));
       Fspring  = max(Fs(springDescription,y),0);
       c = -Fspring + (vLatch.^2).*mLoad./(R.*sqrt(1-(xLatch./R).^2)).*(xLatch.^2./(R.^2.*(1-(xLatch./R).^2))+1)+xLatch.*mLoad./(R.*mLatch.*sqrt(1-(xLatch./R).^2)).*Fl(latchDescription,muk,R,xLatch,y,vLatch,vLoad,[],[]);
       if b == 0
            N = c./(-a);
            aLatch = (Fl(latchDescription,muk,R,xLatch,y,vLatch,vLoad,[],[])+N.*sinTheta-Ff(N,fricInfo).*cosTheta)./mLatch;
       else
            if ~isreal(a)
                a = 0;
            end
            N = fmincon(@(N) findRoot(N, a, b, c, fricInfo), Fspring./cosTheta, [], [], [], [], 0, [], [], options);

            if (Fl(latchDescription,mus,R,xLatch,y,vLatch,vLoad,[],[]) + N.*sinTheta < Ff(N,fricInfo)) && (t == 0 || vLatch < tol)
               error('Latch gets Stuck');
            end
            aLatch = (Fl(latchDescription,muk,R,xLatch,y,vLatch,vLoad,[],[])+N.*sinTheta-Ff(N,fricInfo).*cosTheta)./mLatch;
            if ((abs(aLatch) < tol && vLatch < tol) || (Fl(latchDescription,muk,R,xLatch,y,vLatch,vLoad,[],[]) < tol && vLatch < tol))%|| (vLatch < tol && t ~= 0)

                error('Latch gets stuck');
            end
       end
       dw = [w(2); aLatch];
    else
       dw = [w(2); Fl(latchDescription,muk,R,xLatch,y,vLatch,vLoad,[],[])./mLatch];
    end
end

end

function [diff] = findRoot(N,a,b,c,fricInfo)
diff = ((c+a.*N)./b - Ff(N, fricInfo)).^2;
end
