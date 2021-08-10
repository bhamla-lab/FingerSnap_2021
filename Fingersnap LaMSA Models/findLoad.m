function [FfrMin, FsprMin] = findLoad(theta, fricInfo, initGuess)
options = optimset('Display','off');
zTol = 0.001;
ETol = 10^-3;
maxIts = 25;

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

initGuess = 0;
[N, error, exitflag] = fmincon(@(Nguess) findN(Nguess, theta, fricInfo), initGuess, [], [], [], [], 0, inf, [], options);
FsprMin = N.*cosd(theta)+Ff(N, fricInfo).*sind(theta);
FfrMin = Ff(N, fricInfo);

end

function [err] = findN(Nguess, theta, fricInfo)
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

err = Nguess.*tand(theta)-mus.*Ff(Nguess, fricInfo);
end
