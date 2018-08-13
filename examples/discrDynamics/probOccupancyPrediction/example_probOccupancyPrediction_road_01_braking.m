function completed = example_probOccupancyPrediction_road_01_braking()
% example_probOccupancyPrediction_road_01_braking - example of
% probabilistic occupancy prediction for road traffic when a chain of
% vehicle is predicted and the first one brakes
%
% This example can be found in Sec. 5.5.3 of
% M. Althoff, “Reachability analysis and its application to the safety 
% assessment of autonomous cars”, Dissertation, Technische Universität 
% München, 2010, 
% http://nbn-resolving.de/urn/resolver.pl?urn:nbn:de:bvb:91-diss-20100715-963752-1-4
%
% Syntax:  
%    example_probOccupancyPrediction_road_01_braking
%
% Inputs:
%    no
%
% Outputs:
%    res - boolean 
%
% Example: 
%
% 
% Author:       Matthias Althoff
% Written:      14-October-2009
% Last update:  01-August-2016
% Last revision:---


%------------- BEGIN CODE --------------


%set path
global filePath
filePath = [coraroot '/contDynamics/stateSpaceModels'];
dataPath = [coraroot '/examples/discrDynamics/probOccupancyPrediction/data'];

%load Markov chain of car
cd(dataPath);
file = load('probModel_car_sim_02September2016_newPartitionFormat.mat');
probModel = file.probModel;

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

%load interaction matrix ThetaC
file=load('diss_Theta_opt.mat');
ThetaC=file.ThetaC;

%initial set
initCarA=interval([98; 6], [112; 8]);    %interval for position and speed 

%generate autonomous car trajectory
autonomousCarTrajectory.velocity=3*ones(1,40);
autonomousCarTrajectory.position(1)=120;
for i=1:length(autonomousCarTrajectory.velocity)
    autonomousCarTrajectory.position(i+1)=autonomousCarTrajectory.position(i)...
        +autonomousCarTrajectory.velocity(i)*markovChainSpec.timeStep;
end
autonomousCarTrajectory.mode=ones(1,40)*3;


%set simOptions for braking car
simOptionsA.type='car'; %set as default for simplicity
simOptionsA.runs=16;
simOptionsA.initialStateSet=initCarA;
simOptionsA.stateField=stateField;
simOptionsA.inputField=inputField;
simOptionsA.profileHandle=@profile1;
simOptionsA.transitionMatrix=Tcar;
simOptionsA.interactionMatrix=ThetaC;
simOptionsA.projMat=projMat; 
simOptionsA.GammaFull=GammaFull; 
simOptionsA.mode='autonomousDriving';
simOptionsA.gamma=0.2;
simOptionsA.freeDrivingProb=[0.01 0.04 0.1 0.4 0.4 0.05];
%simOptionsA.initialInputProb=[0 0 0.5 0.5 0 0];
simOptionsA.initialInputProb=[0 0 0 1 0 0];
simOptionsA.autonomousCarTrajectory=autonomousCarTrajectory;
simOptionsA.frontProbVector=[];

	
%instantiate simulation object
carA=simulation(simOptionsA,markovChainSpec);
%carA=simulate(carA); %<--change here
carA=simulateOptimized(carA); %<--change here

%initial set for vehicle B
initCarB=interval([72; 10], [84; 12]); %interval for position and speed
                   
%set simOptions for first follower
simOptionsB=simOptionsA;
simOptionsB.initialStateSet=initCarB;
simOptionsB.frontProbVector=get(carA,'prob');
simOptionsB.mode='vehicleFollowing';

%instantiate simulation object
carB=simulation(simOptionsB,markovChainSpec);
%carB=simulate(carB); %<--change here
carB=simulateOptimized(carB); %<--change here


%initial set for vehicle C (bike)
initCarC=interval([5; 13], [17; 15]); %interval for position and speed 

%change simOptions for third vehicle 
simOptionsC=simOptionsB;
simOptionsC.initialStateSet=initCarC;
simOptionsC.frontProbVector=get(carB,'prob');
simOptionsC.mode='vehicleFollowing';

%instantiate simulation object
carC=simulation(simOptionsC,markovChainSpec);
carC=simulateOptimized(carC); 

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


nrOfSegments=stateField.nrOfSegments;
intervals=stateField.intervals;

%plot position distribution
posField=partition(intervals(1,:),nrOfSegments(1));
for i=[5,10,15]
    figure
    hold on
    plotHisto(posField,posA{i},'k-');
    plotHisto(posField,posB{i},'k--');
    plotHisto(posField,posC{i},'k:');
    xlabel('s [m]');
    ylabel('probability');
end

%plot velocity distribution
velField=partition(intervals(2,:),nrOfSegments(2));
for i=[5,10,15]
    figure
    hold on
    plotHisto(velField,velA{i},'k-');
    plotHisto(velField,velB{i},'k--');
    plotHisto(velField,velC{i},'k:');
    xlabel('v [m/s]');
    ylabel('probability');    
end

u=0; %<-- artificial if statement added 19.04.16 for debugging 
if u==1
%plot input distribution
for i=[5,10,15]
    figure
    hold on
    plotHisto(inputField,inputA{i},'k-');
    plotHisto(inputField,inputB{i},'k--');
    plotHisto(inputField,inputC{i},'k:');
    xlabel('u \in [-1,1]');
    ylabel('probability');    
end
end

%create road
R=road(4,5,1);


%create path
Rstraight=createPath(R,[pi/2,2,0],[0],[40]);

figure;

subplot(1,4,1);
normalizePlot();
%plot velocity distribution
plot(Rstraight,avgVelA,1);
plotCrossing(Rstraight,[36,41]);
axis([-5, 5, 0, 160]);
xlabel('car A');

subplot(1,4,2);
normalizePlot();
%plot velocity distribution
plot(Rstraight,avgVelB,1);
plotCrossing(Rstraight,[36,41]);
axis([-5, 5, 0, 160]);
xlabel('car B');
set(gca,'ytick',[]);

subplot(1,4,3)
normalizePlot();
%plot velocity distribution
plot(Rstraight,avgVelC,1);
plotCrossing(Rstraight,[36,41]);
axis([-5, 5, 0, 160]);
xlabel('car C');
set(gca,'ytick',[]);

subplot(1,4,4)
normalizePlot();



%compute average position probabilities
for iFinal=1:4
    %initialize average position probabilities
    posAvgA{iFinal}=0*posA{1};
    posAvgB{iFinal}=0*posB{1};
    posAvgC{iFinal}=0*posC{1};    
    for i=1:4
        index=4*(iFinal-1)+i;
        posAvgA{iFinal}=posAvgA{iFinal}+posA{index};
        posAvgB{iFinal}=posAvgB{iFinal}+posB{index};
        posAvgC{iFinal}=posAvgC{iFinal}+posC{index};
    end
    %normalize results
    posAvgA{iFinal}=posAvgA{iFinal}/sum(posAvgA{iFinal});
    posAvgB{iFinal}=posAvgB{iFinal}/sum(posAvgB{iFinal});
    posAvgC{iFinal}=posAvgC{iFinal}/sum(posAvgC{iFinal});
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
Rstraight=createPath(dispR,[pi/2,2,0],[0],[40]);

for iStep=1:4
    figure
    %set(gca,'DataAspectRatio',[1 3 3]);    
    
    %plot pA
    subplot(1,3,1);
    plot(Rstraight,posAvgA{iStep},dispDevProbCar,'k');
    %plot crossing
    plotCrossing(Rstraight,[36,41]);
    xlabel('car A');   

    %plot pB
    subplot(1,3,2);
    plot(Rstraight,posAvgB{iStep},dispDevProbCar,'k');
    %plot crossing
    plotCrossing(Rstraight,[36,41]);
    xlabel('car B');
    set(gca,'ytick',[]);    
    
    %plot pC
    subplot(1,3,3);
    plot(Rstraight,posAvgC{iStep},dispDevProbCar,'k');
    %plot crossing
    plotCrossing(Rstraight,[36,41]);  
    xlabel('car C');
    set(gca,'ytick',[]);    
end

%example completed
completed = 1;



function normalizePlot()

%plot lowest and highest value for average probability
%plot using own methods
IH=[0.1 0.2; 0.1 0.2];
V=vertices(IH);
plot(V,'grayTones',0);
plot(V,'grayTones',18);

%------------- END OF CODE --------------
