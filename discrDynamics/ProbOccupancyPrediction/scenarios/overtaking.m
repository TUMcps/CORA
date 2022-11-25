function overtaking()
% changed: 02-October-2009 
% updated: 18-August-2016

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

%load Markov chain of bicycle
[FileName,PathName] = uigetfile('','Select the Markov chain model for the bicycle');
cd(PathName);
file=load(FileName);
probModel_bicycle=file.probModel;

%obtain Tbicycle
Tbicycle = probModel_bicycle.T;

%obtain state and input field
stateField_bicycle = probModel_bicycle.stateField;
inputField_bicycle = probModel_bicycle.inputField;

clc

%initial set
initCarA=interval([2; 14], [8;16]); %interval for position and speed

%generate autonomous car trajectory
autonomousCarTrajectory.velocity=16*ones(1,40);
autonomousCarTrajectory.position(1)=5;
for i=1:length(autonomousCarTrajectory.velocity)
    autonomousCarTrajectory.position(i+1)=autonomousCarTrajectory.position(i)...
        +autonomousCarTrajectory.velocity(i)*markovChainSpec.timeStep;
end
autonomousCarTrajectory.mode=ones(1,40)*3;

%set simOptions
simOptionsA.type='car'; %set as default for simplicity
simOptionsA.runs=18;
simOptionsA.initialStateSet=initCarA;
simOptionsA.stateField=stateField;
simOptionsA.inputField=inputField;
simOptionsA.profileHandle=@profile1;
simOptionsA.transitionMatrix=Tcar;
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

%initial set for vehicle B
initCarB=interval([30; 10], [40;12]); %interval for position and speed
                   
%change simOptions for second vehicle 
simOptionsB=simOptionsA;
simOptionsB.initialStateSet=initCarB;
simOptionsB.profileHandle=@profile3;
simOptionsB.mode='freeDriving';

%instantiate car B
carB=simulation(simOptionsB,markovChainSpec);
carB=simulateOptimized(carB);
posB=get(carB,'posProb');

%change simOptions for second vehicle/ second driving possibility
simOptionsC=simOptionsB;
simOptionsC.profileHandle=@profile4;

%instantiate car C
carC=simulation(simOptionsC,markovChainSpec);
carC=simulateOptimized(carC);
posC=get(carC,'posProb');

%initial set for vehicle C (bike)
initCarD=interval([40; 5], [50; 7]); %interval for position and speed

%change simOptions for third vehicle 
simOptionsD=simOptionsC;
simOptionsD.type='bicycle';
simOptionsD.stateField=stateField_bicycle;
simOptionsD.inputField=inputField_bicycle;
simOptionsD.profileHandle=@profile1;
simOptionsD.initialStateSet=initCarD;
simOptionsD.transitionMatrix=Tbicycle;

%instantiate car D
carD=simulation(simOptionsD,markovChainSpec);
carD=simulateOptimized(carD);
posD=get(carD,'posProb');

% %taking both paths is equally probable:
% for i=1:simOptionsA.runs
%     posB{i}=0.5*posB{i};
%     posC{i}=0.5*posC{i};
% end

%create road
Rcar=road(2,5,4);
Rbicycle=road(3.5,5,7);
dispR=road(4,5,7);

%create path1
R1=createPath(Rcar,[pi/2,2.1,0],[0,0.04,-0.04,0,-0.04,0.04,0],[4,4,4,5,4,4,16]);
dispR1=createPath(dispR,[pi/2,2.1,0],[0,0.04,-0.04,0,-0.04,0.04,0],[4,4,4,5,4,4,16]);
%create path2
R2=createPath(Rcar,[0,-75+8/7,150],[0,-pi/6,0],[12,3,25]);
dispR2=createPath(dispR,[0,-75+8/7,150],[0,-pi/6,0],[12,3,25]);
%create path3
R3=createPath(Rcar,[0,-75+8/7,150],[0,+pi*6/50,+pi/50,0],[12,4,1,24]);
dispR3=createPath(dispR,[0,-75+8/7,150],[0,+pi*6/50,+pi/50,0],[12,4,1,24]);
%create path4
R4=createPath(Rbicycle,[pi/2,2.1,0],[0],[40]);
dispR4=createPath(dispR,[pi/2,2.1,0],[0],[40]);


%deviation probability of the cars
devProbCar=[1,2,2,1];
devProbCar=devProbCar/sum(devProbCar);

deviationFieldCenter=partition([-1, 1],4);
deviationFieldBody=partition([-2, 2],7);
stretch=2; %[m]

%compute body distribution from vehicle center distribution
[dispDevProbCar]=vehicleBodyDistribution(deviationFieldCenter,deviationFieldBody,...
    stretch,devProbCar);


%deviation probability of the bicycle
devProbBicycle=[0,0,0,0,1,2,3];
devProbBicycle=devProbBicycle/sum(devProbBicycle);

deviationFieldCenter=partition([-1.75, 1.75],7);
deviationFieldBody=partition([-2, 2],7);
stretch=0.5; %[m]

%compute body distribution from vehicle center distribution
[dispDevProbBicycle]=vehicleBodyDistribution(deviationFieldCenter,deviationFieldBody,...
    stretch,devProbBicycle);


disp('Intersection time')
%intersection
for iStep=1:length(posA)
    p1(iStep)=intersection(R1,R2,posA{iStep},posB{iStep},...
        devProbCar,devProbCar,probModel.fArray,'car'); 
    p2(iStep)=intersection(R1,R3,posA{iStep},posC{iStep},...
        devProbCar,devProbCar,probModel.fArray,'car');    
    p3(iStep)=intersection(R1,R4,posA{iStep},posD{iStep},...
        devProbCar,devProbBicycle,probModel_bicycle.fArray,'bicycle');
end
p=p1+p2+p3;
t=0:0.5:9;
plot(t,[0,p1],'ok');
hold on
plot(t,[0,p2],'xk');
plot(t,[0,p3],'--k');
plot(t,[0,p]);

%stretch longitudinal probabilities
pathField=partition([0, 200],40);
for iStep=1:length(posA)
    [posA{iStep}]=vehicleBodyDistribution(pathField,pathField,5,posA{iStep});
    [posB{iStep}]=vehicleBodyDistribution(pathField,pathField,5,posB{iStep});
    [posC{iStep}]=vehicleBodyDistribution(pathField,pathField,5,posC{iStep});
    [posD{iStep}]=vehicleBodyDistribution(pathField,pathField,2,posD{iStep});
end


%paths for plotting only
%path5: for straight road
R5=createPath(dispR,[pi/2,2,0],[0],[50]);
RL=createPath(dispR,[0,-15+8/7,148],[-pi/6,0],[3,2],4.165); %lower section
RU=createPath(dispR,[0,-15+8/7,156],[pi/6,0],[3,2],4.165); %upper section
RR=createPath(dispR,[pi/2,4,135],[0],[6]); %right section

%aviobj = avifile('overtakingMovie');
%aviobj=set(aviobj,'fps',3);

for iStep=1:5
    h=figure;
    
    k=4*(iStep-1)+1;
    
%     set(gcf,'Units','normalized');
%     set(h,'position',[0.1,0.1,0.15,0.9]);
%     hold on
%     axis([-5 5 0 245]);
    %set(gca,'DataAspectRatio',[1 3 3]);    
    %plot pA
    plot(dispR1,posA{k},dispDevProbCar,'transB');
    %plot pB 
    plot(dispR2,posB{k},dispDevProbCar,'transR');  
    %plot pC
    plot(dispR3,posC{k},dispDevProbCar,'transR');  
    %plot pD
    plot(dispR4,posD{k},dispDevProbBicycle,'transG');  
    %plot crossing
    plotCrossing(R5,[28,33]);
    plotCrossing(dispR2,[13,40]);
    %plot path of autonomous car
    plotPath(dispR1,1:33); 
    plotPath(RL); 
    plotPath(RU); 
    plotPath(RR); 
    grid
    
    
%    aviobj = addframe(aviobj,getframe);
end

