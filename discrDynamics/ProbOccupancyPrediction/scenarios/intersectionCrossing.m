function intersectionCrossing()
% built: 14-October-2009
% updated: 09-August-2018

%set path
global filePath
filePath = [coraroot '/contDynamics/stateSpaceModels'];

%load Markov chain of car
[FileName,PathName] = uigetfile('','Select the Markov chain model');
cd(PathName);
file=load(FileName);
probModel=file.probModel;
ThetaC=file.ThetaC;

%obtain Tcar, projMat, GammaFull
Tcar = probModel.T;
projMat = probModel.projMat;
GammaFull = probModel.GammaFull;

%Markov-Chain specific settings
markovChainSpec.timeStep = probModel.timeStep;
markovChainSpec.nrOfInputs = prod(probModel.inputField.nrOfSegments);

%obtain state and input field
stateField = probModel.stateField;
inputField = probModel.inputField;

clc

%define initial area with uniform probability distribution
initCarA=interval([40; 10], [52; 12]); %interval for position and velocity
initCarB=interval([5; 10], [17;12]);  %interval for position and velocity
initCarC=interval([25; 6], [37;8]); %interval for position and velocity               
initCarD=interval([5; 8], [17;10]); %interval for position and velocity 
%initCarE=interval([75; 0], [80,1]); %interval for position and velocity 


%set simOptions for braking car
simOptionsA.type='car'; %set as default for simplicity
%simOptionsA.runs=22;
simOptionsA.runs=18;
simOptionsA.initialStateSet=initCarA;
simOptionsA.stateField=stateField;
simOptionsA.inputField=inputField;
simOptionsA.profileHandle=@profile1;
simOptionsA.transitionMatrix=Tcar;
simOptionsA.interactionMatrix=ThetaC;
simOptionsA.projMat=projMat; %<--change here
simOptionsA.GammaFull=GammaFull; %<--change here
simOptionsA.mode='freeDriving';
simOptionsA.gamma=0.2;
simOptionsA.freeDrivingProb=[0.01 0.04 0.1 0.4 0.4 0.05];
%simOptionsA.initialInputProb=[0 0 0.5 0.5 0 0];
simOptionsA.initialInputProb=[0 0 0 1 0 0];
simOptionsA.frontProbVector=[];

	
%instantiate simulation object
carA=simulation(simOptionsA,markovChainSpec);
%carA=simulate(carA); %<--change here
carA=simulateOptimized(carA); %<--change here
posA=get(carA,'posProb');
                   
%set simOptions for first follower
simOptionsB=simOptionsA;
simOptionsB.initialStateSet=initCarB;
simOptionsB.frontProbVector=get(carA,'prob');
simOptionsB.mode='vehicleFollowing';

%instantiate simulation object
carB=simulation(simOptionsB,markovChainSpec);
%carB=simulate(carB); %<--change here
carB=simulateOptimized(carB); %<--change here
posB=get(carB,'posProb');

%compute crossing probability---------------------------
%xSegment=24; %(120 m) ...of other lane
xSegment=20; %(100 m) ...of other lane
intInterval=8; %(4 sec)
%intInterval=6; %(3 sec)
pXD=xIntegration(posA,xSegment,intInterval);
pXC=xIntegration(posB,xSegment,intInterval);
pX=pXC+pXD;
%-------------------------------------------------------

%change simOptions for third vehicle 
simOptionsC=simOptionsB;
simOptionsC.initialStateSet=initCarC;
simOptionsC.frontProbVector=[];
simOptionsC.xSegment=14;
crossProb=max(1-pX,0); %<-- !!!
simOptionsC.tranProb=crossProb;
simOptionsC.mode='roadCrossing';

%instantiate simulation object
carC=simulation(simOptionsC,markovChainSpec);
carC=simulateOptimized(carC); 


%change simOptions for forth vehicle 
simOptionsD=simOptionsC;
simOptionsD.initialStateSet=initCarD;
simOptionsD.frontProbVector=get(carC,'prob');

%instantiate simulation object
carD=simulation(simOptionsD,markovChainSpec);
carD=simulateOptimized(carD); 


%get positio, velocity and input probability distribution
posA=get(carA,'posProb');
velA=get(carA,'velProb');
inputA=get(carA,'inputProb');
avgVelA=get(carA,'avgVel');

posB=get(carB,'posProb');
velB=get(carB,'velProb');
inputB=get(carB,'inputProb');
avgVelB=get(carB,'avgVel');

posC=get(carC,'posProb');
velC=get(carC,'velProb');
inputC=get(carC,'inputProb');
avgVelC=get(carC,'avgVel');

posD=get(carD,'posProb');
velD=get(carD,'velProb');
inputD=get(carD,'inputProb');
avgVelD=get(carD,'avgVel');


nrOfSegments=stateField.nrOfSegments;
intervals=stateField.intervals;

%plot position distribution
posField=partition(intervals(1,:),nrOfSegments(1));
for i=[8,15,simOptionsA.runs]
    figure
    hold on
    plotHisto(posField,posA{i});
    plotHisto(posField,posB{i});
    plotHisto(posField,posC{i});
    plotHisto(posField,posD{i});
    xlabel('s [m]');
    ylabel('probability');
end

%plot velocity distribution
velField=partition(intervals(2,:),nrOfSegments(2));
for i=[8,15,simOptionsA.runs]
    figure
    hold on
    plotHisto(velField,velA{i});
    plotHisto(velField,velB{i});
    plotHisto(velField,velC{i});
    plotHisto(velField,velD{i});
    xlabel('v [m/s]');
    ylabel('probability');    
end

%plot input distribution
for i=[8,15,simOptionsA.runs]
    figure
    hold on
    plotHisto(inputField,inputA{i});
    plotHisto(inputField,inputB{i});
    plotHisto(inputField,inputC{i});
    plotHisto(inputField,inputD{i});
    xlabel('u \in [-1,1]');
    ylabel('probability');    
end


%plot cross probability
t=0.5:markovChainSpec.timeStep:simOptionsA.runs*markovChainSpec.timeStep;
figure
plot(t(1:end),1-crossProb); %<-- !!!


%compute average position probabilities
stepsCombined=4;
for iFinal=1:4
    %initialize average position probabilities
    posAvgA{iFinal}=0*posA{1};
    posAvgB{iFinal}=0*posB{1};
    posAvgC{iFinal}=0*posC{1};   
    posAvgD{iFinal}=0*posD{1};   
    for i=1:stepsCombined
        index=stepsCombined*(iFinal-1)+i;
        posAvgA{iFinal}=posAvgA{iFinal}+posA{index};
        posAvgB{iFinal}=posAvgB{iFinal}+posB{index};
        posAvgC{iFinal}=posAvgC{iFinal}+posC{index};
        posAvgD{iFinal}=posAvgD{iFinal}+posD{index};
    end
    %normalize results
    posAvgA{iFinal}=posAvgA{iFinal}/sum(posAvgA{iFinal});
    posAvgB{iFinal}=posAvgB{iFinal}/sum(posAvgB{iFinal});
    posAvgC{iFinal}=posAvgC{iFinal}/sum(posAvgC{iFinal});
    posAvgD{iFinal}=posAvgD{iFinal}/sum(posAvgD{iFinal});
end


%create road
Rcar=road(2,5,4);
dispR=road(4,5,7);

%deviation probability of the cars
devProbCar=[1,2,2,1];
devProbCar=devProbCar/sum(devProbCar);

deviationFieldCenter=partition([-1, 1],4);
deviationFieldBody=partition([-2, 2],7);
stretch=2; %[m]

%compute body distribution from vehicle center distribution
[dispDevProbCar]=vehicleBodyDistribution(deviationFieldCenter,deviationFieldBody,...
    stretch,devProbCar);

%create path
Rstraight=createPath(dispR,[pi/2,0,5],0,40);
Rcurved=createPath(dispR,[3*pi/2,30,190],[0,-0.1,0.1,0],[12;12;12;4]);
RcurvedLong=createPath(dispR,[3*pi/2,30,190],[0,-0.1,0.1,0],[12;12;12;8]);


    

for iStep=1:4
    figure   
    
    %k=7*(iStep-1)+1;
    k=5*(iStep-1)+1;
    
    %plot pA and pB
    %plot(Rcurved,posAvgA{iStep},dispDevProbCar,'trans');
    %plot(Rcurved,posAvgB{iStep},dispDevProbCar,'trans');
    if iStep==1
        plot(RcurvedLong,posA{k},dispDevProbCar,'trans');
        plot(RcurvedLong,posB{k},dispDevProbCar,'trans');         
    else
        plot(Rcurved,posA{k},dispDevProbCar,'trans');
        plot(Rcurved,posB{k},dispDevProbCar,'trans');    
    end
    %plot crossing
    plotCrossing(Rcurved,[24,27]); 
    
    %plot pC and pD
    %plot(Rstraight,posAvgC{iStep},dispDevProbCar,'trans');
    %plot(Rstraight,posAvgD{iStep},dispDevProbCar,'trans');
    plot(Rstraight,posC{k},dispDevProbCar,'trans');
    plot(Rstraight,posD{k},dispDevProbCar,'trans');    
    %plot crossing 
    plotCrossing(Rstraight,[14,17]); 
    
    %plot remaining lines
    plot([-9.12,-7,-6],[70.36,71,70],'k-');
    plot([-14.09,-10,-6.5,-6],[76.63,80,85,90],'k-');
    plot([4.02,2.2,2],[85.89,85,90],'k-');
    plot([7.65,3,2.2,2],[78.76,75,73,70],'k-');
end


%create road
R=road(4,5,1);

%create path
Rstraight=createPath(R,[pi/2,0,5],0,40);
Rcurved=createPath(R,[3*pi/2,30,190],[0,-0.1,0.1,0],[12;12;12;4]);

%figure;

%subplot(1,5,1);
figure;
normalizePlot();
%plot velocity distribution
plot(Rcurved,avgVelA,1);
plotRoad();
axis([-40, 40, 20, 190]);
% xlabel('car A');

%subplot(1,5,2);
figure;
normalizePlot();
%plot velocity distribution
plot(Rcurved,avgVelB,1);
plotRoad();
axis([-40, 40, 20, 190]);
% xlabel('car B');
% set(gca,'ytick',[]);

%subplot(1,5,3)
figure;
normalizePlot();
%plot velocity distribution
plot(Rstraight,avgVelC,1);
plotRoad();
axis([-40, 40, 20, 190]);
% xlabel('car C');
% set(gca,'ytick',[]);

%subplot(1,5,4)
figure;
normalizePlot();
%plot velocity distribution
plot(Rstraight,avgVelD,1);
plotRoad();
axis([-40, 40, 20, 190]);
% xlabel('car D');
% set(gca,'ytick',[]);


% subplot(1,5,5)
% normalizePlot();
% set(gca,'xtick',[]);
% set(gca,'ytick',[]);


function normalizePlot()

%plot lowest and highest value for average probability
%plot using own methods
IH=[0.1 0.2; 0.1 0.2];
V=vertices(IH);
plot(V,'grayTones',0);
plot(V,'grayTones',18);

function plotRoad()

%create road
R=road(4,5,1);

%create path
Rstraight=createPath(R,[pi/2,0,5],[0],[40]);
Rcurved=createPath(R,[3*pi/2,30,190],[0,-0.1,0.1,0],[12;12;12;4]);

%plot crossing 
plotCrossing(Rcurved,[24,27]); 
plotCrossing(Rstraight,[14,17]); 

%plot remaining lines
plot([-9.12,-7,-6],[70.36,71,70],'k-');
plot([-14.09,-10,-6.5,-6],[76.63,80,85,90],'k-');
plot([4.02,2.2,2],[85.89,85,90],'k-');
plot([7.65,3,2.2,2],[78.76,75,73,70],'k-');
