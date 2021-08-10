function [t,y,dy,d2y,x,dx,d2x,Fspring,netF,unlatchTime,unlatchTimeError,loadUnlatch,timeTO,loadTO,s,Wf,N] = simulation(loadStr,latchStr,startStr,fricStr,detailsStr)
%SIMULATION Takes in information about the sytem to test, returns
%information about the times analyzed; positions, velocities, and
%accelerations of load and latch at each time; spring and net forces at
%each time; unlatch time; error associated with the unlatch time; the
%position and velocity of the load at unlatching; take off time; pposition
%and velocity of the load at take off; arc length travelled; energy lost to
%friction; normal forces between latch and load.
%
%
% Inputs:
%   startStr: structure with info about t =0. fields:
%       startStr.startTheta: double - initial angle between latch and load
%           in degrees
%       startStr.startVx: double - intial horizontal velocity of the latch
% 
%   loadStr: structure with info about spring and loading force. Fields:
%       loadStr.k: double - spring constant
%       loadStr.m: double - mass of the load
% 
%   fricStr: structure with info about friction values/model. Fields:
%       fricStr.form: vector - contains info for the nonlinear friction
%       model [F_f = mu*N(N/N_0)^n] Vector has:
%           n: double - exponential factor N is raised to. Must be -1<n<0
%           N_0: double - N/N_0. Usually 1.
%       fricStr.mus: double - static coefficient of friction.
%       fricStr.muk: double - dynamic coefficient of friction. 
% 
%   latchStr: structure with info about the latch. Fields:
%       latchStr.R: double - radius of the circular latch
%       latchStr.motor: cell - information on unlatching motor. Contains:
%           Motor type: string - "linear_motor", "constant_force"
%           Additional Information
%               If motor type is "linear_motor", cell also needs:
%                   Motor Definitions: vector
%                       max_Force: double - maximum motor force
%                       max_velocity: double - maximum velocity of
%                       unlatching motor
%                       range_of_motion: double - motor range
%                   Braking: boolean
%                       If true, motor stops at end of range of motion
%                       If false, motor continues regardless of range of
%                       motion
% 
%   detailsStr: structure with additional information. Fields:
%       detailsStr.overShoot: double - for each mu, percentage of time 
%       after t_t0 to continue calculating information.  

mus = fricStr.mus;
muk = fricStr.muk;
mLatch = latchStr.m;
mLoad = loadStr.m;
R = latchStr.R;

springDescription = {loadStr.k, loadStr.yEQ};
startX = sin(startStr.startTheta.*pi./180).*latchStr.R;
latchStartConditions = [startX startStr.startVx];

latchDescription = latchStr.motor;
latchDescription{2} = [latchDescription{2} startX];

overShoot = detailsStr.overShoot;
fricFormInfo = fricStr.form;

detailsStr.overShoot = overShoot;

%%Unlatching Phase

[unlatchTime,unlatchTimeError,t1,w1] = find_tUnlatch(R,mus,muk,mLatch,mLoad,latchDescription,springDescription,latchStartConditions,fricFormInfo);
unlatchTime = real(unlatchTime);


timeTakeOff = inf;

options = odeset('RelTol',1e-5,'AbsTol',1e-7);


if unlatchTime ~= 0
    %For successful unlatching
    
    
    %% Phase 1: Unlatching
    t1 = real(t1);
    w1 = real(w1);
    x1 = w1(:,1);
    dx1 = w1(:,2);
    d2x1 = [];
    y1 = [];
    dy1 = [];
    d2y1 = [];
    Fspring1 = []; 
    netF1 = [];
    Wf = 0;
    s = [];
    Ffs = [];
    N1 = [];
    for i = 1:length(t1)
        temp = latchODE(t1(i),w1(i,:),R,mus,muk,mLatch,latchDescription,springDescription,unlatchTime,mLoad,fricFormInfo);
        d2x1 = [d2x1 temp(2)];
        y1 = [y1 yLoad(R,x1(i))];
        dy1 = [dy1 dyLoad(R,x1(i),dx1(i))];
        d2y1 = [d2y1 d2yLoad(R,x1(i),dx1(i),d2x1(i))];
        Fspring1 = [Fspring1 Fs(springDescription,y1(i))];
        netF1 = [netF1 d2y1(i).*mLoad];
        
        sinTheta = x1(i)./(sqrt(x1(i)^2+(R-y1(i)).^2));
        cosTheta = (R-y1(i))./(sqrt(x1(i)^2+(R-y1(i)).^2));
        
        
        s = [s, asin(sinTheta.*pi./180).*R];
        
        a = cosTheta + mLoad.*x1(i).*sinTheta./(mLoad.*R.*sqrt(1-(x1(1)./R).^2));
        b = sinTheta + mLoad.*x1(i).*cosTheta./(mLoad.*R.*sqrt(1-(x1(i)./R).^2));
        c = -Fspring1(i)+(dx1(i)).^2.*mLoad./(R.*sqrt(1-(x1(i)./R).^2)).*(x1(i).^2./(R.^2.*(1-(x1(i)./R).^2))+1)+x1(i).*mLoad.*Fl(latchDescription,muk,R,x1(i),y1(i),dx1(i),dy1(i),[],[])./(R.*mLatch.*sqrt(1-(x1(i)./R).^2));
        
        Ntemp = linspace(0,max(Fspring1(i),0),1000);
        Ff1 = (c+a.*Ntemp)./b;
        fricInfo = [mus, muk, x1(i), dx1(i), y1(i), dy1(i), R, fricFormInfo];
        Ff2 = [];
        for m = 1:length(Ntemp)
            Ff2 = [Ff2, Ff(Ntemp(m), fricInfo)];
        end
        [~, index] = min(Ff2-Ff1);
        Na = Ntemp(index);
        
        left = netF1(i) - Fspring1(i);
        right = [];
        for m = 1:length(Ntemp)
            right = [right -(Ntemp(m).*cosTheta + Ff(Ntemp(m), fricInfo).*sinTheta)];
        end
        [~, index2] = min(abs(left-right));
        Nb = Ntemp(index2);
        
        N1 = [N1 Nb];
        Ffs = [Ffs Ff(Ntemp(index), fricInfo)];
    
    end
    
    mask = imag(y1) ~= 0;
    t1(mask) = [];
    y1(mask) = [];
    dy1(mask) = [];
    d2y1(mask) = [];
    x1(mask) = [];
    dx1(mask) = [];
    d2x1(mask) = [];
    Fspring1(mask) = [];
    netF1(mask) = [];
    N1(mask) = [];
    Ffs(mask) = [];
    s(mask) = [];
    
    loadUnlatch = [y1(end),dy1(end)];
    latchUnlatch = [x1(end),dx1(end)];
    
    if length(s) > 1
        Wf = cumtrapz((s-s(1)),Ffs);
    else
        Wf = [0];
    end
    
    
    if Fs(springDescription,y1(end)) >= 0
        %If unlatching position finishes after the spring equilibrium
        %position
        
        
        %% Phase 2: Standard harmonic motion
        k = springDescription{1};
        yEQ = springDescription{2};
        wn = sqrt(k./mLoad);
        A = sqrt((loadUnlatch(1)-yEQ).^2+(loadUnlatch(2)./wn).^2);
        phi = atan(loadUnlatch(2)./(wn.*(-loadUnlatch(1)+yEQ)));
        phi = phi.*pi./180;

        latchStartConditions = [x1(end), dx1(end)];    
        timeTO = atan(wn.*(yEQ-loadUnlatch(1))/loadUnlatch(2))/wn+unlatchTime;

        t2 = linspace(unlatchTime,timeTO,100);
        x2 = [];
        dx2 = [];
        d2x2 = [];
        y2 = [];
        dy2 = [];
        d2y2 = [];
        Fspring2 = []; 
        netF2 = [];
        N2 = [];

        for i = 1:length(t2)
            calcY = (loadUnlatch(1)-yEQ).*cos(wn.*(t2(i)-unlatchTime)) + loadUnlatch(2)./wn.*sin(wn.*(t2(i)-unlatchTime))+yEQ;
            y2 = [y2, calcY];
            calcdY = -wn.*(loadUnlatch(1)-yEQ).*sin(wn.*(t2(i)-unlatchTime)) + loadUnlatch(2).*cos(wn.*(t2(i)-unlatchTime));
            dy2 = [dy2 calcdY];
            calcd2Y = -wn.^2.*(loadUnlatch(1)-yEQ).*cos(wn.*(t2(i)-unlatchTime)) - loadUnlatch(2).*wn.*sin(wn.*(t2(i)-unlatchTime));
            d2y2 = [d2y2 calcd2Y];

            if i == 1
                lF = Fl(latchDescription,mus,R,x1(end),y1(end),dx1(end),dy1(end),dx1(end),dy1(end));
            else
                lF = Fl(latchDescription,mus,R,x2(i-1),y2(i-1),dx2(i-1),dy2(i-1),dx2(i-1),dy2(i-1));
            end
            d2x2 = [d2x2 lF./mLoad];
            dx2 = [dx2 lF./mLoad.*(t2(i)-unlatchTime)+latchUnlatch(2)];
            x2 = [x2 lF./(mLoad.*2).*(t2(i)-unlatchTime).^2+latchUnlatch(2).*(t2(i)-unlatchTime)+latchUnlatch(1)];

            Fspring2 = [Fspring2 Fs(springDescription,y2(i))];
            netF2 = [netF2 d2y2(i).*mLoad];        
            N2 = [N2, 0];
        end

        t2 = t2';
        x2 = x2';
        dx2 = dx2';
        d2x2 = d2x2;
        y2 = y2;
        dy2 = dy2;
        d2y2 = d2y2;
        Fspring2 = Fspring2; 
        netF2 = netF2;


        loadTO = [y2(end),dy2(end)];
        latchStartConditions = [x2(end), dx2(end)];

        %% Phase 3: Load mass detatched from spring
        
        t3 = linspace(timeTO, timeTO*overShoot);
        x3 = [];
        dx3 = [];
        d2x3 = [];
        y3 = [];
        dy3 = [];
        d2y3 = [];
        Fspring3 = []; 
        netF3 = [];
        N3 = [];

        for i = 1:length(t3)
            if i == 1
                lF = Fl(latchDescription,mus,R,x2(end),y2(end),dx2(end),dy2(end),dx2(end),dy2(end));
            else
                lF = Fl(latchDescription,mus,R,x3(i-1),y3(i-1),dx3(i-1),dy3(i-1),dx2(i-1),dy3(i-1));
            end

            d2x3 = [d2x3; lF./mLoad];
            dx3 = [dx3; lF./mLoad.*(t3(i)-unlatchTime)+latchUnlatch(2)];
            x3 = [x3 lF./(mLoad.*2).*(t3(i)-unlatchTime).^2+latchUnlatch(2).*(t3(i)-unlatchTime)+latchUnlatch(1)];

            d2y3 = [d2y3 0];
            dy3 = [dy3 loadTO(2)];
            y3 = [y3 loadTO(2).*(t3(i)-timeTO)+loadTO(1)];
            Fspring3 = [Fspring3 Fs(springDescription,y3(i))];
            netF3 = [netF3 d2y3(i).*mLoad];
            N3 = [N3 0];
        end
        t3 = t3';
        x3 = x3';
        dx3 = dx3;
        d2x3 = d2x3';

        loadUnlatch = [y1(end),dy1(end),d2y1(end)];
        t = [t1; t2; t3];
        x = [x1; x2; x3];
        dx = [dx1; dx2; dx3];
        d2x = [d2x1'; d2x2'; d2x3'];
        y = [y1'; y2'; y3'];
        dy = [dy1'; dy2'; dy3'];
        d2y = [d2y1'; d2y2'; d2y3'];
        Fspring = [Fspring1'; Fspring2'; Fspring3']; 
        netF = [netF1'; netF2'; netF3'];
        N = [N1'; N2'; N3'];
    else
        loadUnlatch = [y1(end),dy1(end),d2y1(end)];
        t = t1;
        x = x1;
        dx = dx1;
        d2x = d2x1';
        y = y1';
        dy = dy1';
        d2y = d2y1';
        Fspring = Fspring1';
        netF = netF1';
        N = N1';
        timeTO = unlatchTime;
        loadTO = loadUnlatch(1:2);
    end
    
    
else
    % If unlatching is unsuccessful, or if the latch is a constant velocity
    % latch
    unlatchTime = NaN;
    [w1,rest] = strtok('latchDescription');
    if strcmpi(w1,'constant_velocity')
        unlatchTime = inf;
        unlatchTimeError = 0;
        loadUnlatch = [NaN,NaN];
        timeTO = inf;
        loadTO = [NaN,NaN];

        [w1, rest] = strtok(latchDescription);

        temp = str2num(rest)./2*Rx;


        if temp ~= 0
            [t,w] = ode45(@(t,w) latchODE(t,w,R,mus,muk,mLatch,latchDescription,springDescription,unlatchTime,timeTakeOff), [0 temp], latchStartConditions);
            x = w(:,1);
            dx = w(:,2);
            d2x = NaN;
            y = [];
            dy = [];
            d2y = [];
            Fspring = []; 
            netF = [];

            for i = 1:length(t)
                temp = latchODE(t(i),w(i,:),R,mus,muk,mLatch,latchDescription,springDescription,unlatchTime,timeTakeOff);
                d2x = [d2x temp(2)];
                y = [y yLoad(R,t(i),x(i),unlatchTime,timeTakeOff)];
                dy = [dy dyLoad(R,t(i),x(i),dx(i),unlatchTime,timeTakeOff)];
                d2y = [d2y d2yLoad(R,t(i),x(i),dx(i),d2x(i),unlatchTime,timeTakeOff)];
                Fspring = [Fspring Fs(springDescription,y(i))];
                netF = [netF d2y(i).*mLoad];
            end
            d2x = d2x';
            y = y';
            dy = dy';
            d2y = d2y';
            Fspring = Fspring';
            netF = netF';
            s = 0;
            Wf = 0;
            N= [];
        else
            [t,w] = ode45(@(t,w) latchODE(t,w,Rx,Ry,n,mus,muk,mLatch,latchDescription,springDescription,unlatchTime,timeTakeOff), [0 temp+overShoot], latchStartConditions);
            x = w(:,1);
            dx = w(:,2);
            d2x = [];
            y = [];
            dy = [];
            d2y = [];
            Fspring = []; 
            netF = [];

            for i = 1:length(t)
                temp = latchODE(t(i),w(i,:),Rx,Ry,n,mus,muk,mLatch,latchDescription,springDescription,unlatchTime,timeTakeOff);
                d2x = [d2x temp(2)];
                y = [y yLoad(R,t(i),x(i),unlatchTime,timeTakeOff)];
                dy = [dy dyLoad(R,t(i),x(i),dx(i),unlatchTime,timeTakeOff)];
                d2y = [d2y d2yLoad(R,t(i),x(i),dx(i),d2x(i),unlatchTime,timeTakeOff)];
                Fspring = [Fspring Fs(springDescription,y(i))];
                netF = [netF d2y(i).*mLoad];
            end
            d2x = d2x';
            y = y';
            dy = dy';
            d2y = d2y';
            Fspring = Fspring';
            netF = netF';
            s = 0;
            Wf = 0;
            N = [];
            
        end
    else
        xTemp = latchStartConditions(1);
        yTemp = yLoad(R,xTemp);
        
        t = [];
        x = [];
        dx = [];
        d2x = [];
        y = [];
        dy = [];
        d2y = NaN;
        Fspring = []; 
        netF = [];
        unlatchTimeError = NaN;
        loadUnlatch = [NaN,NaN, NaN];
        timeTO = NaN;
        loadTO = [NaN,NaN];
        s = NaN;
        Wf = NaN;
        N = [];
        
        timeTO = unlatchTime;
    end

end



end

function [diff] = findRoot(N,a,b,c,fricInfo)
diff = ((c+a.*N)./b - Ff(N, fricInfo)).^2;

end