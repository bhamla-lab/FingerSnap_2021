clear
%close all
%%
%The purpose of this case is to explore the effects of friction on a
%circular latch being loaded by the phenomenological model of force
%loading.

%Latch Parameters
R = 0.1;

%Load Parameters
k = 1000;
    %Spring constant
mLoad = 1;
    %Mass of the load
    
muRatio = 1;
    %Ratio of static to dynamic friction
    
%Misc Parameters

startTheta0 = 0; %degrees
startV = 0;
    %Describes initial position and velocity of latch 
    
mLatch = 1;
    %latch mass in kg
unlatchMotor = {'linear_motor' [100,0.05,0.1] false};
    %First value in the cell describes the motor type.
    %Second value includes details on the motor. For a linear motor:
        %First value in vector is the maximum force of the unlatching motor
        %Second value in vector is the maximum motor unlatching velocity
        %Third value in the vector is the motor range of motion
    %Third value indicates whether the latch continues motion after the
    %maximum range of motion is reached - termed "braking"

fricformInfo = [0, 1];
%First number represents nonlinearity n (should be between 0 and 1), second
%number represents N_0 should usually be 1).

alpha = 250;
    %Factor determining how much the loaded force increases by with respect
    %to the coefficient of friction
F0 = 40;
    %Initial spring force applied when coefficient of friction is 0

loading = 'phenomenological';

mu_s = [linspace(0,1,100)];

overShoot = 10;
    %Describes percentage of time to overshoot take off time by.
       
loadStr.k = k;
loadStr.loading = loading;
loadStr.loadingInfo = {alpha, F0};
loadStr.loadMaxF = 100;
loadStr.m = mLoad;

startStr.startTheta = startTheta0;
startStr.startVx = 0;

fricStr.form = fricformInfo;
fricStr.mus = mu_s;
fricStr.mu_rat = muRatio; 
    
latchStr.R = R;
latchStr.motor = unlatchMotor;
latchStr.m = mLatch;

detailsStr.overShoot = overShoot;
detailsStr.backVec = [13, 26, 53, 56, 60, 67];

[kineticsa] = fricVar(startStr,loadStr,fricStr,latchStr,detailsStr);