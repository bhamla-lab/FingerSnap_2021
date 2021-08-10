function [Ff] = Ff(N,fricInfo)
%FF Calculates horizontal component of friction on latch from coefficient
%of friction, description of spring, Ry, normal force, and latch and load
%positions. 
%   First calls function N to find normal force at latch position x and
%   load position y. Then multiplies value by cosTheta to find component of
%   N perpendicular to surface. Then multiplies this value by mu to
%   find the force of Friction between the surfaces. Finally, multiplies
%   this by cosTheta to find the horizontal component of friction that acts
%   on the latch. Equations for sinTheta and cosTheta are based on latch
%   geometry.
mus = fricInfo(1);
muk = fricInfo(2);
xLatch = fricInfo(3);
vLatch = fricInfo(4);
yLoad = fricInfo(5);
vLoad = fricInfo(6);
R = fricInfo(7);
if length(fricInfo) > 7
    otherfricInfo = fricInfo(8:end);
end

sinTheta = xLatch./(sqrt(xLatch^2+(R-yLoad).^2));
cosTheta = (R-yLoad)./(sqrt(xLatch^2+(R-yLoad).^2)); 

if vLoad == 0
    Ff = mus*N*(N/otherfricInfo(2)).^otherfricInfo(1);
else
    Ff = muk*N*(N/otherfricInfo(2)).^otherfricInfo(1);
end

end

