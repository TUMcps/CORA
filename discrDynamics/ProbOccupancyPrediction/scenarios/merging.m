function merging()
% changed: 02-October-2009 
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

%initial set
initCarA=interval([3;10], [3;10]); %interval for the position and velocity

%generate autonomous car trajectory
autonomousCarTrajectory.velocity=[10,8,ones(1,8)*6,7,8,9,10,11,12,13,14];
autonomousCarTrajectory.position(1)=0;
for i=1:length(autonomousCarTrajectory.velocity)
    autonomousCarTrajectory.position(i+1)=autonomousCarTrajectory.position(i)...
        +autonomousCarTrajectory.velocity(i)*markovChainSpec.timeStep;
end
autonomousCarTrajectory.mode=[2,2,3*ones(1,8),4,4,4,4,4,4,4,4]; 

%set simOptions
simOptionsA.type='car'; %set as default for simplicity
simOptionsA.runs=15;
simOptionsA.initialStateSet=initCarA;
simOptionsA.stateField=stateField;
simOptionsA.inputField=inputField;
simOptionsA.profileHandle=@profile1;
simOptionsA.transitionMatrix=Tcar;
simOptionsA.interactionMatrix=ThetaC;
simOptionsA.projMat=projMat; %<--change here
simOptionsA.GammaFull=GammaFull; %<--change here
simOptionsA.mode='autonomousDriving';
simOptionsA.gamma=0.2;
simOptionsA.freeDrivingProb=[0.01 0.04 0.1 0.4 0.4 0.05];
simOptionsA.initialInputProb=[0 0 0 1 0 0];
simOptionsA.autonomousCarTrajectory=autonomousCarTrajectory;
simOptionsA.frontProbVector=[];

	
%instantiate simulation object
carA=simulation(simOptionsA,markovChainSpec);
carA=simulateOptimized(carA);
posA=get(carA,'posProb');
pTotalA=get(carA,'total');

%initial set for vehicle B
initCarB=interval([10;14], [20;15]); %interval for the position and velocity
                   
%change simOptions for second vehicle 
simOptionsB=simOptionsA;
simOptionsB.initialStateSet=initCarB;
simOptionsB.mode='freeDriving';

%instantiate car B
carB=simulation(simOptionsB,markovChainSpec);
carB=simulateOptimized(carB);
posB=get(carB,'posProb');
pTotalB=get(carB,'total');

%initial set for vehicle C
initCarC=interval([80;5], [90;7]);  %interval for the position and velocity

%change simOptions for third vehicle
simOptionsC=simOptionsB;
simOptionsC.initialStateSet=initCarC;
%simOptionsC.initialInputProb=[0 0 1 0 0 0];

%instantiate car C
carC=simulation(simOptionsC,markovChainSpec);
carC=simulateOptimized(carC);
posC=get(carC,'posProb');
pTotalC=get(carC,'total');

%initial set for vehicle C (bike)
initCarD=interval([20; 12], [30;13]);  %interval for the position and velocity

%change simOptions for forth vehicle 
simOptionsD=simOptionsC;
simOptionsD.initialStateSet=initCarD;
simOptionsD.mode='vehicleFollowing';
simOptionsD.frontProbVector=get(carC,'prob');


%instantiate car D
carD=simulation(simOptionsD,markovChainSpec);
carD=simulateOptimized(carD);
posD=get(carD,'posProb');
pTotalD=get(carD,'total');

%intersection comes afterwards...

%create road
R=road(2,5,4);
dispR=road(4,5,7);

%create path1
R1=createPath(R,[0,-25+8/7+0.3,30],[0,+pi*6/50,+pi/50,0],[2,4,1,23]);
dispR1=createPath(dispR,[0,-25+8/7+0.3,30],[0,+pi*6/50,+pi/50,0],[2,4,1,23]);
%create path2
R2=createPath(R,[-pi/2,-2,140],[0],[40]);
dispR2=createPath(dispR,[-pi/2,-2,140],[0],[40]);
%create path3
R3=createPath(R,[pi/2,2,0],[0],[40]);
dispR3=createPath(dispR,[pi/2,2,0],[0],[40]);

%deviation probability of the cars
devProbCar=[1,2,2,1];
devProbCar=devProbCar/sum(devProbCar);

deviationFieldCenter=partition([-1, 1],4);
deviationFieldBody=partition([-2, 2],7);
stretch=2; %[m]

%compute body distribution from vehicle center distribution
[dispDevProbCar]=vehicleBodyDistribution(deviationFieldCenter,deviationFieldBody,...
    stretch,devProbCar);

disp('Intersection geometric')

%intersection
for iStep=1:length(posA)
    %old intersection method
    p1(iStep)=intersection(R1,R2,posA{iStep},posB{iStep},...
        devProbCar,devProbCar,probModel.fArray,'car'); 
    p2(iStep)=intersection(R1,R3,posA{iStep},posC{iStep},...
        devProbCar,devProbCar,probModel.fArray,'car');    
    p3(iStep)=intersection(R1,R3,posA{iStep},posD{iStep},...
        devProbCar,devProbCar,probModel.fArray,'car');
end

disp('Intersection dynamic')
for iStep=1:length(posA)    
    %intersection method by Alexander
    pA.pos=posA{iStep};
    pA.total=pTotalA.T{iStep};
    pB.pos=posB{iStep};
    pB.total=pTotalB.T{iStep};
    pC.pos=posC{iStep};
    pC.total=pTotalC.T{iStep};
    pD.pos=posD{iStep};
    pD.total=pTotalD.T{iStep};
end
%compute crash probabilities
p=p1+p2+p3;

%plot crash probabilities
t=0:0.5:(simOptionsA.runs*0.5);
plot(t,[0,p],'k-');

%stretch longitudinal probabilities
pathField=partition([0, 200],40);
for iStep=1:length(posA)
    [posA{iStep}]=vehicleBodyDistribution(pathField,pathField,5,posA{iStep});
    [posB{iStep}]=vehicleBodyDistribution(pathField,pathField,5,posB{iStep});
    [posC{iStep}]=vehicleBodyDistribution(pathField,pathField,5,posC{iStep});
    [posD{iStep}]=vehicleBodyDistribution(pathField,pathField,5,posD{iStep});
end


%paths for plotting only
%path4: for straight road
R4=createPath(dispR,[pi/2,2,0],[0],[40]);
R5=createPath(dispR,[0,-25+8/7,30],[0],[3]);
RL=createPath(dispR,[0,-15+8/7,28],[-pi/6,0],[3,2],4.165); %lower section
RU=createPath(dispR,[0,-15+8/7,36],[pi/6,0],[3,2],4.165); %upper section
RR=createPath(dispR,[pi/2,4,20],[0],[6]); %right section

% aviobj = avifile('mergingMovie');
% aviobj=set(aviobj,'fps',3);

for iStep=1:length(posA)
    figure
    
%     set(gcf,'Units','normalized');
%     set(h,'position',[0.1,0.1,0.15,0.9]);
%     hold on
%     axis([-5 5 0 245]);
    %set(gca,'DataAspectRatio',[1 3 3]);    
    %plot pA
    plot(dispR1,posA{iStep},dispDevProbCar,'trans');
    %plot pB
    plot(dispR2,posB{iStep},dispDevProbCar,'trans');
    %plot pC
    plot(dispR3,posC{iStep},dispDevProbCar,'trans'); 
    %plot pD
    plot(dispR3,posD{iStep},dispDevProbCar,'trans');
    %plot crossing
    plotCrossing(R4,[5,9]);
    plotCrossing(R5,[3,3]);
    %plot path of autonomous car
    plotPath(R1,1:29);  

    
    plotPath(RL);
    plotPath(RU);
    plotPath(RR);
    grid
    changeAxis();
%     aviobj = addframe(aviobj,getframe);
end

function changeAxis()

axis([-10 5 0 150]);
xlabel('');
ylabel('');

