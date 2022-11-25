function knownBehavior()
% built: 22-October-2009
% update: 09-August-2018


%Markov-Chain specific settings
markovChainSpec.timeStep=0.5;
markovChainSpec.nrOfInputs=6;

%initialize vehicle to get state and input field
[~,~,stateField]=initCar(3);


%generate autonomous car trajectory
autonomousCarTrajectory.velocity=[5 6 7 8 9 10 11 11.5 12 12 11.5 10 8 7 7 7 8 9 9 9 9];
autonomousCarTrajectory.position(1)=10;
for i=1:length(autonomousCarTrajectory.velocity)
    autonomousCarTrajectory.position(i+1)=autonomousCarTrajectory.position(i)...
        +autonomousCarTrajectory.velocity(i)*markovChainSpec.timeStep;
end
autonomousCarTrajectory.mode=ones(1,40)*3;



%get field
field=stateField;

%position, velocity and mode vector
posVec=autonomousCarTrajectory.position;
velVec=autonomousCarTrajectory.velocity;

% map results to uniform probabilities
% assume that car is x meters long
% assume that speed is uncertain within y m/s
uncertainPos=3;
uncertainVel=1;

iStep=10 %<-- set iStep
    
%generate interval hull for time point solution
IHtp=interval([posVec(iStep)-uncertainPos;velVec(iStep)-uncertainVel],...
    [posVec(iStep)+uncertainPos;velVec(iStep)+uncertainVel]);

%get interval hull for time interval
minPos=min(posVec(iStep),posVec(iStep+1));
maxPos=max(posVec(iStep),posVec(iStep+1));

minVel=min(velVec(iStep),velVec(iStep+1));
maxVel=max(velVec(iStep),velVec(iStep+1));    

IHti=interval([minPos-uncertainPos; minVel-uncertainVel],...
    [maxPos+uncertainPos; maxVel+uncertainVel]);

%convert intervals to probabilities
%time point
[~, pTP] = exactIntersectingCells(field,IHtp);
%time interval
[~, pTI] = exactIntersectingCells(field,IHti);


%plots---------------------------------------------------------------------
figure
%plot state field
plot(field);

%plot trajectory
plot(posVec(1:20),velVec(1:20),'ko-');

%plot uncertain region
%plot(IHti,[1 2],'blackEdge');
plot(IHti,[1 2],'k-','lineWidth',2);
%set labels
xlabel('s [m]');
ylabel('v [m/s]');
%set axis limits
axis([0, 100, 2, 16]);

%plot distribution in new field
figure
plot(field);
%dummy Markov chain forplotting
MC = markovchain(field);
plotP(MC,pTI,'k');
%set labels
xlabel('s [m]');
ylabel('v [m/s]');
%set axis limits
axis([0, 100, 2, 16]);

