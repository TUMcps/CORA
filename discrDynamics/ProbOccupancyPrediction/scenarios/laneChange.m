function laneChange()
% updated: 10-August-2018

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
initCarA=interval([5;14], [15;16]); %interval for position and velocity
initCarB=interval([72;14], [84;16]);... %interval for position and velocity
initCarC=interval([72;4], [84;8]); %interval for position and velocity 
initCarD=interval([30;14], [40;16]); %interval for position and velocity
                   

%generate autonomous car trajectory
autonomousCarTrajectory.velocity=15*ones(1,40);
autonomousCarTrajectory.position(1)=10;
for i=1:length(autonomousCarTrajectory.velocity)
    autonomousCarTrajectory.position(i+1)=autonomousCarTrajectory.position(i)...
        +autonomousCarTrajectory.velocity(i)*markovChainSpec.timeStep;
end
autonomousCarTrajectory.mode=ones(1,40)*3;

%set simOptions
simOptionsA.type='car'; %set as default for simplicity
simOptionsA.runs=10;
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
%simOptionsA.initialInputProb=[0 0 0 1 0 0];
simOptionsA.initialInputProb=[0 0 0.5 0.5 0 0];
simOptionsA.autonomousCarTrajectory=autonomousCarTrajectory;
simOptionsA.frontProbVector=[]; 

	
%instantiate simulation object
carA=simulation(simOptionsA,markovChainSpec);
carA=simulateOptimized(carA);

                   
%change simOptions for second vehicle 
simOptionsB=simOptionsA;
simOptionsB.initialStateSet=initCarB;
simOptionsB.mode='freeDriving';

%instantiate car B
carB=simulation(simOptionsB,markovChainSpec);
carB=simulateOptimized(carB);



%change simOptions for third vehicle 
simOptionsC=simOptionsB;
simOptionsC.initialStateSet=initCarC;
simOptionsC.initialInputProb=[1 0 0 0 0 0];
simOptionsC.mode='freeDriving';

%instantiate car C
carC=simulation(simOptionsC,markovChainSpec);
carC=simulateOptimized(carC);



%change simOptions for forth vehicle 
simOptionsD=simOptionsC;
simOptionsD.initialStateSet=initCarD;
simOptionsD.frontProbVector{1}=get(carC,'prob');
simOptionsD.frontProbVector{2}=get(carB,'prob');
simOptionsD.frontProbVector{3}=get(carA,'prob');
simOptionsD.initialInputProb=[0 0 0 1 0 0];
simOptionsD.mode='laneChanging';

%instantiate car D
carD=simulation(simOptionsD,markovChainSpec);
carD=simulateOptimized(carD);
lcEvolProb=get(carD,'lcEvolProb');


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



%create road
R=road(4,5,1);

%create path
Rleft=createPath(R,[pi/2,-2,0],[0],[40]);
Rright=createPath(R,[pi/2,2,0],[0],[40]);

figure;

subplot(1,5,1);
normalizePlot();
%plot velocity distribution
plot(Rright,avgVelA,1);
plotRoad;
xlabel('car A');

subplot(1,5,2);
normalizePlot();
%plot velocity distribution
plot(Rright,avgVelB,1);
plotRoad;
xlabel('car B');
set(gca,'ytick',[]);

subplot(1,5,3)
normalizePlot();
%plot velocity distribution
plot(Rleft,avgVelC,1);
plotRoad;
xlabel('car C');
set(gca,'ytick',[]);

subplot(1,5,4)
normalizePlot();
%plot velocity distribution
plot(Rleft,avgVelD.left,1);
plot(Rright,avgVelD.right,1);
plotRoad;
xlabel('car D');
set(gca,'ytick',[]);

subplot(1,5,5)
normalizePlot();


%create road
R=road(4,5,8); %<--outer segments have to be chosen as 0
dispR=road(4,5,7);

%create path
Rleft=createPath(R,[pi/2,-2,0],[0],[40]);
Rright=createPath(R,[pi/2,2,0],[0],[40]);

dispRleft=createPath(dispR,[pi/2,-2,0],[0],[40]);
dispRright=createPath(dispR,[pi/2,2,0],[0],[40]);


%deviation probability of the vehicles
[devProbA,dispDevProbA]=deviationProbability(Rright,length(posA));
[devProbB,dispDevProbB]=deviationProbability(Rright,length(posA));
[devProbC,dispDevProbC]=deviationProbability(Rleft,length(posA));
[devProbDleft,dispDevProbDleft,devProbDright,dispDevProbDright]=deviationProbability(Rleft,lcEvolProb);


for iStep=1:5
    figure;
    hold on  
    
    k=2*(iStep-1)+1;
    
    %plot pA, pB and pC 
    plot(dispRright,posA{k},dispDevProbA(k,:),'trans');
    plot(dispRright,posB{k},dispDevProbB(k,:),'trans');
    plot(dispRleft,posC{k},dispDevProbC(k,:),'trans');
    
    %plot pD
    %get maximum probability for normalization
    pMaxLeft=max(posD.left{k})*1.1*max(dispDevProbDleft(k,:));
    pMaxRight=max(posD.right{k})*1.1*max(dispDevProbDright(k,:));
    pMax=max(pMaxLeft,pMaxRight);
    normVal=(sum(posD.left{k})+sum(posD.right{k}))/pMax;    
    plot(dispRleft,posD.left{k},dispDevProbDleft(k,:),'trans',normVal);   
    plot(dispRright,posD.right{k},dispDevProbDright(k,:),'trans',normVal);    
      
    %plot road
    plotRoad;
  
end


% for iStep=1:length(posA)
%     figure;
%     hold on
%     %set(gca,'DataAspectRatio',[1 3 3]);    
%     %plot pA
%     plot(dispRright,posA{iStep},dispDevProbA(iStep,:),'b');
%     %plot pB
%     plot(dispRright,posB{iStep},dispDevProbB(iStep,:),'g');    
%     %plot pC
%     plot(dispRleft,posC{iStep},dispDevProbC(iStep,:),'r'); 
%     %plot pD
%     %get maximum probability for normalization
%     pMaxLeft=max(posD.left{iStep})*1.1*max(dispDevProbDleft(iStep,:));
%     pMaxRight=max(posD.right{iStep})*1.1*max(dispDevProbDright(iStep,:));
%     pMax=max(pMaxLeft,pMaxRight);
%     normVal=(sum(posD.left{iStep})+sum(posD.right{iStep}))/pMax;    
%     plot(dispRleft,posD.left{iStep},dispDevProbDleft(iStep,:),'k',normVal);   
%     plot(dispRright,posD.right{iStep},dispDevProbDright(iStep,:),'k',normVal); 
%     %plot crossing
%     plotCrossing(Rright,[1,1]);
% end



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
Rleft=createPath(R,[pi/2,-2,0],[0],[40]);
Rright=createPath(R,[pi/2,2,0],[0],[40]);

%plot crossing 
plotCrossing(Rright,[40,41]) 

