%clear
%close all
%%
%The purpose of this case is to explore the effects of latch, load, and spring paramaters on the resulting
%kinematics of the system while the circular latch is moving with a constant force being applied to it. 

%Latch Parameters
R = 0.1;


%Load Parameters
k = 1000;
yEQ = 0.01;
springDescription = {k,yEQ}; 
    %First number represents k in N/m, second number represents initial displacement (yEq)
mLoad = 1; 
    %Mass of the load in kg
    
    
muRatio = 1;
    
%Misc Parameters

startTheta0 = 0;
latchStartConditions = [sin(startTheta0.*pi./180).*R 0];
    %Describes initial position and velocity of latch in m and m/s.
    %If latch is moving at a constant velocity, then second value (initial velocity) has to be equal to latch
    %velocity
    
mLatch = 1; 
    %latch mass in kg
latchKinematic = 1;
latchDescription = {'linear_motor' [100,0.05,0.1,latchStartConditions(1)] false};
%latchDescription = {'linear_motor' [90,2000,20] false};
    %If constant force, latchDescription = 'constant_force F' where F is force acting on Latch in Newtons
    
    %If constant velocity, latchDescription = 'constant_velocity v' where v is velocity of Latch in m/s
    
    
overShoot = 10;
    %Describes percentage of time to overshoot take off time by.
    
    


bar = jet;


numPoints1 = 75;

fricformInfo = [0, 1];
tMat = zeros(numPoints1,1);
vMat = zeros(numPoints1,1);
sMat = zeros(numPoints1,1);
vTOMat = zeros(numPoints1,1);
WfMat = zeros(numPoints1,1);
a0Mat = zeros(numPoints1,1);
Fsprings = zeros(numPoints1,1);
Ff0s = zeros(numPoints1,1);
aMaxMat = zeros(numPoints1,1);

mu_ss = [];

minmu_s = 0.01;
maxmu_s =1;
yEQ0 = yEQ;

alpha = 250;
F0 = 40;

loading = 'phenomenological';

FsprMinMax = 100;

%mu_s = [linspace(0,0.12406,100), linspace(0.12407,0.150217,3), linspace(0.150218,0.197064,100), linspace(0.24453,1,70)];
mu_s = [linspace(0,0.3,50), linspace(0.301,1,25)];
%mu_s = logspace(log10(minmu_s),log10(maxmu_s),numPoints1);
mu_k = mu_s./muRatio;
mu_ss = [mu_ss; mu_s];


kinetics = [];
kineticsa = [];
firstReached = false;
reached = true;
maxReached = false;
for i = 1:numPoints1  
    fricInfo = [mu_s(i), mu_s(i)./1.5, R-R*cosd(startTheta0), 0, latchStartConditions(1), 0, R, fricformInfo];
    
    switch loading
        case 'phenomenological'
            FsprMin= alpha.*mu_s(i)+F0;
        case 'nonlinear'
            [FfrMin, FsprMin] = findLoad(startTheta0, fricInfo, 10^8);
            FsprMin = FsprMin;
    end
    if FsprMin > FsprMinMax
        FsprMin = FsprMinMax;
    end

    startX = sin(startTheta0.*pi./180).*R;
    latchStartConditions = [startX 0];
    
    yEQ = FsprMin/k + (R-R.*cosd(startTheta0));
    springDescription = {k, yEQ};

    tic
    [t,yLoads,vLoads,aLoads,~,~,~,~,netF,unlatchTime,unlatchTimeError,loadUnlatch,timeTO,loadTO,s,Wf,N]...
        = simulation(R,mu_k(i),mu_s(i),mLatch,mLoad,latchDescription,springDescription,latchStartConditions,overShoot,fricformInfo);

    toc
    tMat(i) = unlatchTime;
    vMat(i) = loadUnlatch(2);
    sMat(i) = loadUnlatch(1);
    vTOMat(i) = loadTO(2);
    WfMat(i) = Wf(end);
    a0Mat(i) = aLoads(1); 
    Fsprings(i) = FsprMin;
    %Ff0s(i,m) = FfrMin;
    aMaxMat(i) = loadUnlatch(3); 
    
    if isnan(timeTO)
        maxReached = true;
    end
    
    kineticsa = [kineticsa struct('mu',mu_s(i),'t',t,'y',yLoads,'dy',vLoads,'N',N,'s',s,'Wf',Wf,'netF',netF,'tul',unlatchTime,'tTO',timeTO)];
    if i==60 %~firstReached && FsprMin~=FsprMinMax && FsprMin~=0
        mu1 = struct('mu',mu_s(i),'t',t,'y',yLoads,'dy',vLoads,'N',N,'tTO',timeTO);
        if ~isnan(unlatchTime)
            firstReached = true;
        end
        
    elseif FsprMin==FsprMinMax && reached
        mu2 = struct('mu',mu_s(i),'t',t,'y',yLoads,'dy',vLoads,'N',N,'tTO',timeTO);
        reached = false;
    end
    if ~maxReached
        mu3 = struct('mu',mu_s(i),'t',t,'y',yLoads,'dy',vLoads,'N',N,'tTO',timeTO);
    end
    
end

% ogmu_ss = mu_ss;
% kinetics = [mu1 mu2 mu3];
% ogtMat = tMat;
% ogvMat = vMat;
% ogsMat = sMat;
% ogvTOMat = sMat;
% %vTOMat = smooth(vTOMat);
% ogWfMat = WfMat;
% oga0Mat = a0Mat; 
% ogFsprings= Fsprings;
% %Ff0s(i,m) = FfrMin;
% ogaMaxMat = aMaxMat;
% ogkineticsa = kineticsa;
%%
% mu_ss = ogmu_ss;
% kinetics = [mu1 mu2 mu3];
% tMat = ogtMat;
% vMat = ogvMat;
% sMat = ogsMat;
% vTOMat = ogsMat;
% %vTOMat = smooth(vTOMat);
% WfMat = ogWfMat;
% a0Mat = oga0Mat; 
% Fsprings= ogFsprings;
% %Ff0s(i,m) = FfrMin;
% aMaxMat = ogaMaxMat;
% kineticsa = ogkineticsa;



set(0,'defaulttextinterpreter','latex')
fSize = 16;

cF = pwd;
folderName = sprintf('Density- (%d), Loading- (%s,%d), Linear Motor- (maxF=%d,maxV=%.2d,range=%.2d), Spring- (k=%d), Latch- (m=%d,R=%.2d), Friction- (%.2d,%.2d), Load- (m=%.2d), Theta- (%d)',...
    numPoints1, loading(1:2), FsprMinMax, latchDescription{2}(1),latchDescription{2}(2),latchDescription{2}(3), k, mLatch, R, fricformInfo(1),fricformInfo(2),mLoad, startTheta0);

mkdir(folderName)
old = cd(folderName);
%save('workspace');

[~, maxVidx]  = max(vTOMat);
idx1 = 1;
idx2 = maxVidx;
idx3 = numPoints1-5;

kinFig = figure('Name','Kinetics','WindowState', 'maximized');
subplot(3,1,1)
hold on
plot(real(kineticsa(idx1).t),real(kineticsa(idx1).y),'DisplayName',['\mu=' num2str(kineticsa(idx1).mu)])
plot(real(kineticsa(idx2).t),real(kineticsa(idx2).y),'DisplayName',['\mu=' num2str(kineticsa(idx2).mu)])
plot(real(kineticsa(idx3).t),real(kineticsa(idx3).y),'DisplayName',['\mu=' num2str(kineticsa(idx3).mu)])
xlim([min(real(kineticsa(idx1).t)) 1.8*max([kineticsa(idx1).tTO kineticsa(idx2).tTO kineticsa(idx3).tTO])])
%ylim([0, max([real(kineticsa(idx1).y(end)) real(kineticsa(idx2).y(end)) real(kineticsa(idx3).y(end))])])
set(gca, 'XTickLabel', [])
legend
ylabel('$y_{m} (m)$','FontSize', fSize)
hold off

subplot(3,1,2)
hold on
plot(real(kineticsa(idx1).t),real(kineticsa(idx1).dy))
plot(real(kineticsa(idx2).t),real(kineticsa(idx2).dy))
plot(real(kineticsa(idx3).t),real(kineticsa(idx3).dy))
xlim([min(real(kineticsa(idx1).t)) 1.8*max([kineticsa(idx1).tTO kineticsa(idx2).tTO kineticsa(idx3).tTO])])
set(gca, 'XTickLabel', [])
ylabel('$v_{m} (m/s)$','FontSize', fSize)
hold off

subplot(3,1,3)
hold on
plot(real(kineticsa(idx1).t),real(kineticsa(idx1).N))
plot(real(kineticsa(idx2).t),real(kineticsa(idx2).N))
plot(real(kineticsa(idx3).t),real(kineticsa(idx3).N))
xlim([min(real(kineticsa(idx1).t)) 1.8*max([kineticsa(idx1).tTO kineticsa(idx2).tTO kineticsa(idx3).tTO])])
xlabel('$t (s)$','FontSize', fSize)
ylabel('$N (N)$','FontSize', fSize)
hold off
saveas(kinFig,'Kinetics')


eidx1 = 1;
eidx2 = idx3;
E1Fig = figure('Name','Energy Distribution','WindowState', 'maximized');
subplot(1,2,1)
hold on
U1 = [];
for i = 1:length(real(kineticsa(eidx1).y))
    if real(kineticsa(eidx1).y(i)) < Fsprings(eidx1)./k
        U1 = [U1; 0.5.*k.*(Fsprings(eidx1)./k-real(kineticsa(eidx1).y(i))).^2];
    else
        U1 = [U1;0];
    end
    
end
[~, unlatchIndex] = min(abs(kineticsa(eidx1).t-kineticsa(eidx1).tul));
unlatchIndex = unlatchIndex+2;
[~, toIndex] = min(abs(kineticsa(eidx1).t-kineticsa(eidx1).tTO));
K1 = 0.5.*mLoad.*(real(kineticsa(eidx1).dy)).^2;
Eloss1 = U1(1)-U1-K1;
Estruc1 = Eloss1(1:length(kineticsa(eidx1).Wf))-kineticsa(eidx1).Wf;

area(kineticsa(eidx1).t,Eloss1, 'FaceColor', 'm', 'EdgeAlpha', 0,  'FaceAlpha', 0.5)
area(kineticsa(eidx1).t,Eloss1'-[kineticsa(eidx1).Wf,zeros(1,length(kineticsa(eidx1).t)-length(kineticsa(eidx1).Wf))+kineticsa(eidx1).Wf(end)],'FaceColor', 'c', 'EdgeAlpha', 0,  'FaceAlpha', 0.5)
plot(real(kineticsa(eidx1).t),U1,'r','DisplayName','U')
plot(real(kineticsa(eidx1).t),K1,'k','DisplayName','K')
plot(real(kineticsa(eidx1).t),Eloss1,'b','DisplayName','E_d')
ylabel('$E (J)$','FontSize', fSize)
xlabel('$t (s)$','FontSize', fSize)
hold off

subplot(1,2,2)
hold on
U2 = [];
for i = 1:length(real(kineticsa(eidx2).y))
    if real(kineticsa(eidx2).y(i)) < Fsprings(eidx2)./k
        U2 = [U2; 0.5.*k.*(Fsprings(eidx2)./k-real(kineticsa(eidx2).y(i))).^2];
    else
        U2 = [U2;0];
    end
    
end
[~, unlatchIndex] = min(abs(kineticsa(eidx2).t-kineticsa(idx1).tul));
unlatchIndex = unlatchIndex+2;
[~, toIndex] = min(abs(kineticsa(eidx2).t-kineticsa(eidx2).tTO));
K2 = 0.5.*mLoad.*(real(kineticsa(eidx2).dy)).^2;
Eloss2 = U2(1)-U2-K2;
Estruc2 = Eloss2(1:length(kineticsa(eidx2).Wf))-kineticsa(eidx2).Wf;

area(kineticsa(eidx2).t,Eloss2, 'FaceColor', 'm', 'EdgeAlpha', 0,  'FaceAlpha', 0.5)
area(kineticsa(eidx2).t,Eloss2'-[kineticsa(eidx2).Wf,zeros(1,length(kineticsa(eidx2).t)-length(kineticsa(eidx2).Wf))+kineticsa(eidx2).Wf(end)],'FaceColor', 'c', 'EdgeAlpha', 0,  'FaceAlpha', 0.5)
plot(real(kineticsa(eidx2).t),U2,'r','DisplayName','U')
plot(real(kineticsa(eidx2).t),K2,'k','DisplayName','K')
plot(real(kineticsa(eidx2).t),Eloss2,'b','DisplayName','E_d')

hold off
saveas(E1Fig,'Energetic Evolution')


mask = isnan(tMat) ~= 0 | imag(tMat) ~= 0;
mu_ss(mask) = [];
tMat(mask) = [];
vMat(mask) = [];
sMat(mask) = [];
vTOMat(mask) = [];
%vTOMat = smooth(vTOMat);
WfMat(mask) = [];
a0Mat(mask) = []; 
Fsprings(mask) = [];
%Ff0s(i,m) = FfrMin;
aMaxMat(mask) = [];
kineticsa(mask) = [];

% toremove = [51:55];
% mu_ss(toremove) = [];
% tMat(toremove) = [];
% vMat(toremove) = [];
% sMat(toremove) = [];
% vTOMat(toremove) = [];
% %vTOMat = smooth(vTOMat);
% WfMat(toremove) = [];
% a0Mat(toremove) = []; 
% Fsprings(toremove) = [];
% %Ff0s(i,m) = FfrMin;
% aMaxMat(toremove) = [];
% kineticsa(toremove) = [];


trendFig = figure('Name','Friction trends','WindowState', 'maximized');
title('Linear Motor (90, 0.05, 0.1, false), Spring(k=1000,a=220,F0=40), latch(m=1,R=0.1),friction(-0.75,1), load(m=1)')
subplot(2,2,1)
plot(mu_ss,Fsprings,'k-','LineWidth',2)
ylim([0,max(Fsprings)])
xlim([0, max(mu_ss)])
ylabel('$F_{s,0} (N)$','FontSize', fSize)
set(gca, 'XTickLabel', [])

subplot(2,2,3)
plot(mu_ss,real([smooth(tMat(1:maxVidx),1); smooth(tMat((maxVidx+1):end),1)]),'k-','LineWidth',2)
ylim([0,max(real(tMat))])
xlim([0, max(mu_ss)])
ylabel('$t_{ul} (s)$','FontSize', fSize)
xlabel('$\mu$','FontSize', fSize)



smoothedVTO = [real(smooth(vTOMat(1:maxVidx),5)); real(smooth(vTOMat((maxVidx+1):end),7))];

subplot(2,2,[2, 4])
plot(mu_ss,smoothedVTO,'k-','LineWidth',2)
ylim([0,max(real(vTOMat))])
xlim([0, max(mu_ss)])
ylabel('$v_{TO} (m/s)$','FontSize', fSize)
xlabel('$\mu$','FontSize', fSize)

saveas(trendFig,'Trends')


EFig = figure('Name','Energetics','WindowState', 'maximized');

KE = 0.5.*mLoad.*real((smoothedVTO)).^2;
%KE = 0.5.*mLoad.*real(smooth(vTOMat)).^2;
PE = 0.5.*Fsprings.^2./k;
Ed = smooth(PE-KE,3);

subplot(2,1,1)
hold on
plot(mu_ss,KE,'k-','LineWidth',3)
plot(mu_ss,PE,'r-','LineWidth',2)
plot(mu_ss,Ed,'b-','LineWidth',2)
xlim([0, max(mu_ss)])
ylabel('$E (J)$','FontSize', fSize)
set(gca, 'XTickLabel', [])
hold off

subplot(2,1,2)
dPE = gradient(PE,mu_ss);
dEd = [smooth(gradient(Ed(1:maxVidx),mu_ss(1:maxVidx)),3); smooth(gradient(Ed((maxVidx+1):end),mu_ss((maxVidx+1):end)),9)];
hold on
plot(mu_ss,dPE,'r-','LineWidth',2)
xlim([0, max(mu_ss)])
plot(mu_ss,dEd,'b-','LineWidth',2)
hold off
xlabel('$\mu$','FontSize', fSize)
ylabel('$\frac{dE}{d\mu}$','FontSize', fSize)

ax2 = axes('Position', [0.4, 0.2, 0.2, 0.2]);
hold on
plot(mu_ss,dPE,'r-','LineWidth',2)
plot(mu_ss,dEd,'b-','LineWidth',2)
xlim([mu_ss(maxVidx-1) mu_ss(maxVidx+2)])
ylim([min([min(dPE((maxVidx-1):(maxVidx+2))) min(dEd((maxVidx-1):(maxVidx+2)))]) max([max(dPE((maxVidx-1):(maxVidx+2))) max(dEd((maxVidx-1):(maxVidx+2)))])])
hold off

ax3 = axes('Position', [0.7, 0.2, 0.2, 0.2]);
hold on
diff = dPE-dEd;
plot(mu_ss,diff,'k-','LineWidth',2)
plot(mu_ss,zeros(1,length(mu_ss)),'-','Color',[0.25 0.25 0.25])
xlim([mu_ss(maxVidx-2) mu_ss(maxVidx+2)])
ylim([min(diff((maxVidx-1):(maxVidx+2))) max(diff((maxVidx-1):(maxVidx+2)))])
ylabel('$\frac{dU_{sp,0}}{d\mu}-\frac{dE_{d,tot}}{d\mu}$','FontSize', fSize-2)
hold off

saveas(EFig, 'Energetics')

cd(old)

%saveas(f, 'Linear Motor (90, 0.05, 0.1, false), Spring(k=1000,a=220,F0=40), latch(m=1,R=0.1), friction(-0.75,1), load(m=1)')
