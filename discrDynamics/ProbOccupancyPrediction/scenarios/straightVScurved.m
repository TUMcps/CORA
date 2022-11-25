function straightVScurved()
% built: 09-October-2009
% updated: 09-August-2016


%set path
global filePath
filePath = [coraroot '/contDynamics/stateSpaceModels'];

%load Markov chain of car
[FileName,PathName] = uigetfile('','Select the Markov chain model');
cd(PathName);
file=load(FileName);
probModel=file.probModel;

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

%initial set
initCarA=interval([2; 12], [8;14]);  %interval for position and velocity

%set simOptions for straight road
simOptionsA.type='car'; %set as default for simplicity
simOptionsA.runs=20;
simOptionsA.initialStateSet=initCarA;
simOptionsA.stateField=stateField;
simOptionsA.inputField=inputField;
simOptionsA.profileHandle=@profile1;
simOptionsA.transitionMatrix=Tcar;
simOptionsA.projMat=projMat; %<--change here
simOptionsA.GammaFull=GammaFull; %<--change here
simOptionsA.mode='freeDriving';
simOptionsA.gamma=0.2;
simOptionsA.freeDrivingProb=[0.01 0.04 0.1 0.4 0.4 0.05];
%simOptionsA.initialInputProb=[0 0 0.5 0.5 0 0];
simOptionsA.initialInputProb=[0 0 0 1 0 0];
simOptionsA.frontProbVectorSeries=[];

%set simOptions for curved road
simOptionsB=simOptionsA;
simOptionsB.profileHandle=@profile4;
	
%instantiate simulation object
carA=simulation(simOptionsA,markovChainSpec);
%carA=simulate(carA); %<--change here
carA=simulateOptimized(carA); %<--change here

%instantiate simulation object
carB=simulation(simOptionsB,markovChainSpec);
profile on
%carB=simulate(carB); %<--change here
carB=simulateOptimized(carB); %<--change here
profile off
profile viewer

%get positio, velocity and input probability distribution
posA=get(carA,'posProb');
velA=get(carA,'velProb');
inputA=get(carA,'inputProb');
avgVelA=get(carA,'avgVel');

posB=get(carB,'posProb');
velB=get(carB,'velProb');
inputB=get(carB,'inputProb');
avgVelB=get(carB,'avgVel');


nrOfSegments=stateField.nrOfSegments;
intervals=stateField.intervals;

%plot position distribution
posField=partition(intervals(1,:),nrOfSegments(1));
for i=16
    figure
    hold on
    plotHisto(posField,posA{i},'k-');
    plotHisto(posField,posB{i},'k--');
    xlabel('s [m]');
    ylabel('probability');
end

%plot velocity distribution
velField=partition(intervals(2,:),nrOfSegments(2));
for i=16
    figure
    hold on
    plotHisto(velField,velA{i},'k-');
    plotHisto(velField,velB{i},'k--');
    xlabel('v [m/s]');
    ylabel('probability');  
end

%plot input distribution
for i=16
    figure
    hold on
    plotHisto(inputField,inputA{i},'k-');
    plotHisto(inputField,inputB{i},'k--');
    xlabel('u \in [-1,1]');
    ylabel('probability');      
end


%create road
R=road(4,5,1);


%straight road-------------------------------------------------------------

%create path
Rstraight=createPath(R,[pi/2,2,0],[0],[80]);

figure

%plot lowest and highest value for average probability
%plot using own methods
IH=[0.1 0.2; 0.1 0.2];
V=vertices(IH);
plot(V,'grayTones',4);
plot(V,'grayTones',18);

%plot velocity distribution
plot(Rstraight,avgVelA,1);

%plot crossing
plotCrossing(Rstraight,[43,80]);
axis([-10, 10, 0, 200]);


%curved road---------------------------------------------------------------

%create path
Rcurved=createPath(R,[pi/2,2,0],[0,+pi*6/50,+pi/50,0],[12,4,1,24]);

figure

%plot lowest and highest value for average probability
%plot using own methods
IH=[0.1 0.2; 0.1 0.2];
V=vertices(IH);
plot(V,'grayTones',4);
plot(V,'grayTones',18);

%plot velocity distribution
plot(Rcurved,avgVelB,1);

%plot crossing
plotCrossing(Rcurved,[35,80]);
axis([-60, 10, 0, 90]);
%--------------------------------------------------------------------------

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
Rstraight=createPath(dispR,[pi/2,2,0],[0],[80]);
Rcurved=createPath(dispR,[pi/2,2,0],[0,+pi*6/50,+pi/50,0],[12,4,1,24]);

for iStep=16
    figure
    %set(gca,'DataAspectRatio',[1 3 3]);    
    %plot pA
    plot(Rstraight,posA{iStep},dispDevProbCar,'k');
    %plot crossing
    plotCrossing(Rstraight,[43,80]);
    axis([-10, 10, 0, 200]);
    
    figure
    %plot pB
    plot(Rcurved,posB{iStep},dispDevProbCar,'k');
    %plot crossing
    plotCrossing(Rcurved,[35,80]); 
    axis([-60, 10, 0, 90]);
end

