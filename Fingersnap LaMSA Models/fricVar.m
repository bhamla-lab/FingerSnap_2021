function [kineticsa] = fricVar(startStr,loadStr,fricStr,latchStr,detailsStr)
%FRICVAR Takes in information about the sytem as well as values of mu_s to
%test, creates a folder with a txt document describing the parameters
%tested, creates and saves four folders, and outputs a structure with
%fields containing the experiment results.
%Calls simulation for each mu_s value specified in the fricStr structure
%variable input.
%
%Inputs:
%   startStr: structure with info about t =0. fields:
%       startStr.startTheta: double - initial angle between latch and load
%           in degrees
%       startStr.startVx: double - intial horizontal velocity of the latch
%
%   loadStr: structure with info about spring and loading force. Fields:
%       loadStr.k: double - spring constant
%       loadStr.loading: string - type of loading - "nonlinear" or
%       "phenomenological"
%       loadStr.loadingInfo: cell - additional information
%           If loading is phenomenological, cell requires:
%               alpha: double - linear relationship between loaded force
%               and mu
%               F0: double - loaded force at mu=0
%           If loading is nonlinear, cell can be empty.
%       loadStr.loadMaxF: double - maximum force spring can be loaded to 
%       loadStr.m: double - mass of the load
%
%   fricStr: structure with info about friction values/model. Fields:
%       fricStr.form: vector - contains info for the nonlinear friction
%       model [F_f = mu*N(N/.N_0)^n] Vector has:
%           n: double - exponential factor N is raised to. Must be -1<n<0
%           N_0: double - N/N_0. Usually 1.
%       fricStr.mus: vector - contains all mu_s which should be tested.
%       fricStr.mu_rat: double - ratio between mu_s and mu_k. 
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
%       detailsStr.backVec: double - indexes of mu_vec to graph on the
%       'Kinetics' figure.       

detailsStr.numPoints = length(fricStr.mus);

startX = sin(startStr.startTheta.*pi./180).*latchStr.R;
latchStartConditions = [startX startStr.startVx];

startY = yLoad(latchStr.R,startX);
startvY = dyLoad(latchStr.R,startX, startStr.startVx);


tMat = zeros(detailsStr.numPoints,1);
vMat = zeros(detailsStr.numPoints,1);
sMat = zeros(detailsStr.numPoints,1);
vTOMat = zeros(detailsStr.numPoints,1);
WfMat = zeros(detailsStr.numPoints,1);
a0Mat = zeros(detailsStr.numPoints,1);
Fsprings = zeros(detailsStr.numPoints,1);
aMaxMat = zeros(detailsStr.numPoints,1);

mu_ss = [];

maxLoadF = loadStr.loadMaxF;
mu_ss = [mu_ss; fricStr.mus];

kineticsa = [];
reached = true;
maxReached = false;
for i = 1:detailsStr.numPoints
    tic
    muk = fricStr.mus(i)./fricStr.mu_rat;
    mus = fricStr.mus(i);
    fricInfo = [mus, muk, startX, startStr.startVx, startY, startvY, latchStr.R, fricStr.form];
    
    switch loadStr.loading
        case 'phenomenological'
            loadFspr= loadStr.loadingInfo{1}.*mus+loadStr.loadingInfo{2};
        case 'nonlinear'
            [FfrMin, loadFspr] = findLoad(startStr.startTheta, fricInfo, 10^8);
    end
    if loadFspr > maxLoadF
        loadFspr = maxLoadF;
    end
    
    yEQ = loadFspr/loadStr.k + (latchStr.R-latchStr.R.*cosd(startStr.startTheta));
    springDescription = {loadStr.k, yEQ};
    loadStr.yEQ = yEQ;
    tempFricStr.mus = mus;
    tempFricStr.muk = muk;
    tempFricStr.form = fricStr.form;
    
    
    [t,yLoads,vLoads,aLoads,~,~,~,~,netF,unlatchTime,unlatchTimeError,loadUnlatch,timeTO,loadTO,s,Wf,N]...
        = simulation(loadStr,latchStr,startStr,tempFricStr,detailsStr);

    
    tMat(i) = unlatchTime;
    vMat(i) = loadUnlatch(2);
    sMat(i) = loadUnlatch(1);
    vTOMat(i) = loadTO(2);
    WfMat(i) = Wf(end);
    a0Mat(i) = aLoads(1); 
    Fsprings(i) = loadFspr;
    aMaxMat(i) = loadUnlatch(3); 
    
    if isnan(timeTO)
        maxReached = true;
    end
    
    kineticsa = [kineticsa struct('mu',mus,'t',t,'y',yLoads,'dy',vLoads,'N',N,'s',s,'Wf',Wf,'netF',netF,'tul',unlatchTime,'tTO',timeTO)];
    if i==60 
        mu1 = struct('mu',mus,'t',t,'y',yLoads,'dy',vLoads,'N',N,'tTO',timeTO);
        if ~isnan(unlatchTime)
            firstReached = true;
        end
        
    elseif loadFspr==maxLoadF && reached
        mu2 = struct('mu',mus,'t',t,'y',yLoads,'dy',vLoads,'N',N,'tTO',timeTO);
        reached = false;
    end
    if ~maxReached
        mu3 = struct('mu',mus,'t',t,'y',yLoads,'dy',vLoads,'N',N,'tTO',timeTO);
    end
    pointTime = toc;
    fprintf("Point %d/%d completed in %.2d s.\n",i,detailsStr.numPoints,pointTime)
end

%%
files = dir;
dirFlags = [files.isdir];
subFolders = files(dirFlags);
maxFile = 0; 
for g = 1:length(subFolders)
    %if subFolders(g).name
        maxFile = max([maxFile, str2double(subFolders(g).name)]);
end
maxFile = maxFile;

set(0,'defaulttextinterpreter','latex')
fSize = 16;
minG = 1;
maxG = 225;
cF = pwd;
folderName = num2str(maxFile+1);

mkdir(folderName)
old = cd(folderName);

fid = fopen('Experiment Properties.txt','wt');
fprintf(fid,'Point Density: %d\nLoading: Loading Type=%s, Maximum Loading Motor Force= %d\nUnlatching Motor: Maximum Motor Force=%d, Maximum Motor Velocity=%.2d, Motor Range=%.2d\nSpring: k=%d\nLatch: m=%d, R=%.2d\nFriction: n = %.2d, mu range: %.2d to %.2d\nLoad: m=%.2d,\nTheta_0: %d degrees\n',...
    detailsStr.numPoints, loadStr.loading, maxLoadF, latchStr.motor{2}(1),latchStr.motor{2}(2),latchStr.motor{2}(3), loadStr.k, latchStr.m, latchStr.R, fricStr.form(1),min(fricStr.mus),max(fricStr.mus),loadStr.m,startStr.startTheta);
fclose(fid);

[~, maxVidx]  = max(vTOMat);
idx1 = 1;
idx2 = maxVidx;
idx3 = detailsStr.numPoints-1;
backroundVectors = detailsStr.backVec;


kinFig = figure('Name','Kinetics','WindowState', 'maximized');
subplot(3,1,1)
hold on
for i = backroundVectors(end:-1:1)
    curG =  maxG-(maxG-minG)./(fricStr.mus(end)-fricStr.mus(1)).*fricStr.mus(i);
    plot(real(kineticsa(i).t),real(kineticsa(i).y),'Color',[curG,curG,curG]./255,'LineWidth',1,'DisplayName',['\mu=' num2str(fricStr.mus(i))])
end
plot(real(kineticsa(idx1).t),real(kineticsa(idx1).y),'DisplayName',['\mu=' num2str(kineticsa(idx1).mu)])
plot(real(kineticsa(idx2).t),real(kineticsa(idx2).y),'DisplayName',['\mu=' num2str(kineticsa(idx2).mu)])
plot(real(kineticsa(idx3).t),real(kineticsa(idx3).y),'DisplayName',['\mu=' num2str(kineticsa(idx3).mu)])
xlim([min(real(kineticsa(idx1).t)) 1.8*max([kineticsa(idx1).tTO kineticsa(idx2).tTO kineticsa(idx3).tTO])])
set(gca, 'XTickLabel', [])
legend
ylabel('$y_{m} (m)$','FontSize', fSize)
hold off

subplot(3,1,2)
hold on
for i = backroundVectors(end:-1:1)
    curG =  maxG-(maxG-minG)./(fricStr.mus(end)-fricStr.mus(1)).*fricStr.mus(i);
    plot(real(kineticsa(i).t),real(kineticsa(i).dy),'Color',[curG,curG,curG]./255,'LineWidth',1,'DisplayName',['\mu=' num2str(fricStr.mus(i))])
end
plot(real(kineticsa(idx1).t),real(kineticsa(idx1).dy))
plot(real(kineticsa(idx2).t),real(kineticsa(idx2).dy))
plot(real(kineticsa(idx3).t),real(kineticsa(idx3).dy))
xlim([min(real(kineticsa(idx1).t)) 1.8*max([kineticsa(idx1).tTO kineticsa(idx2).tTO kineticsa(idx3).tTO])])
set(gca, 'XTickLabel', [])
ylabel('$v_{m} (m/s)$','FontSize', fSize)
hold off

subplot(3,1,3)
hold on
for i = backroundVectors(end:-1:1)
    curG =  maxG-(maxG-minG)./(fricStr.mus(end)-fricStr.mus(1)).*fricStr.mus(i);
    plot(real(kineticsa(i).t),real(kineticsa(i).N),'Color',[curG,curG,curG]./255,'LineWidth',1,'DisplayName',['\mu=' num2str(fricStr.mus(i))])
end
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
    if real(kineticsa(eidx1).y(i)) < Fsprings(eidx1)./loadStr.k
        U1 = [U1; 0.5.*loadStr.k.*(Fsprings(eidx1)./loadStr.k-real(kineticsa(eidx1).y(i))).^2];
    else
        U1 = [U1;0];
    end
    
end
[~, unlatchIndex] = min(abs(kineticsa(eidx1).t-kineticsa(eidx1).tul));
unlatchIndex = unlatchIndex+2;
[~, toIndex] = min(abs(kineticsa(eidx1).t-kineticsa(eidx1).tTO));
K1 = 0.5.*loadStr.m.*(real(kineticsa(eidx1).dy)).^2;
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
    if real(kineticsa(eidx2).y(i)) < Fsprings(eidx2)./loadStr.k
        U2 = [U2; 0.5.*loadStr.k.*(Fsprings(eidx2)./loadStr.k-real(kineticsa(eidx2).y(i))).^2];
    else
        U2 = [U2;0];
    end
    
end
[~, unlatchIndex] = min(abs(kineticsa(eidx2).t-kineticsa(idx1).tul));
unlatchIndex = unlatchIndex+2;
[~, toIndex] = min(abs(kineticsa(eidx2).t-kineticsa(eidx2).tTO));
K2 = 0.5.*loadStr.m.*(real(kineticsa(eidx2).dy)).^2;
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
WfMat(mask) = [];
a0Mat(mask) = []; 
Fsprings(mask) = [];
aMaxMat(mask) = [];
kineticsa(mask) = [];

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



smoothedVTO = [real(smooth(vTOMat(1:maxVidx),1)); real(smooth(vTOMat((maxVidx+1):end),1))];

subplot(2,2,[2, 4])
plot(mu_ss,smoothedVTO,'k-','LineWidth',2)
ylim([0,max(real(vTOMat))])
xlim([0, max(mu_ss)])
ylabel('$v_{TO} (m/s)$','FontSize', fSize)
xlabel('$\mu$','FontSize', fSize)

saveas(trendFig,'Trends')


EFig = figure('Name','Energetics','WindowState', 'maximized');

KE = 0.5.*loadStr.m.*real((smoothedVTO)).^2;
PE = 0.5.*Fsprings.^2./loadStr.k;
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
dEd = [smooth(gradient(Ed(1:maxVidx),mu_ss(1:maxVidx)),1); smooth(gradient(Ed((maxVidx+1):end),mu_ss((maxVidx+1):end)),1)];
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
xlim([mu_ss(maxVidx-1) mu_ss(maxVidx+1)])
ylim([min([min(dPE((maxVidx-1):(maxVidx+2))) min(dEd((maxVidx-1):(maxVidx+2)))]) max([max(dPE((maxVidx-1):(maxVidx+2))) max(dEd((maxVidx-1):(maxVidx+2)))])])
hold off

ax3 = axes('Position', [0.7, 0.2, 0.2, 0.2]);
hold on
diff = dPE-dEd;
plot(mu_ss,diff,'k-','LineWidth',2)
plot(mu_ss,zeros(1,length(mu_ss)),'-','Color',[0.25 0.25 0.25])
xlim([mu_ss(maxVidx-1) mu_ss(maxVidx+1)])
ylim([min(diff((maxVidx-1):(maxVidx+2))) max(diff((maxVidx-1):(maxVidx+2)))])
ylabel('$\frac{dU_{sp,0}}{d\mu}-\frac{dE_{d,tot}}{d\mu}$','FontSize', fSize-2)
hold off

saveas(EFig, 'Energetics')

cd(old)


end

