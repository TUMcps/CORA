function res = testLongDuration_markovchain_carReach
% testLongDuration_markovchain_carReach - unit_test_function for building a
% Markov chain of a car using simulation results
%
% Syntax:  
%    res = testLongDuration_markovchain_carReach
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
% Written:      31-July-2016
% Last update:  ---
% Last revision:---


%------------- BEGIN CODE --------------

%set path
global filePath
filePath = [coraroot '/contDynamics/stateSpaceModels'];

%load car model
[HA,options,stateField,inputField,changeSpeed] = testAuxiliaryFct_initCar;

%set gamma value (how often inputs change)
gamma = 0.2;

%set final time and time steps
finalTime=0.5;
timeSteps=5;

%Initialize Markov Chain
MC=markovchain(stateField);

%obtain number of segments of the state discretization
nrOfSegments=get(stateField,'nrOfSegments');

%total number of discrete inputs, states, positions and velocities
totalNrOfInputs=prod(get(inputField,'nrOfSegments'));
totalNrOfPositions=nrOfSegments(1);
totalNrOfVelocities=nrOfSegments(2);

% check for even number of inputs
if rem(totalNrOfInputs,2)>0
        disp('Number of Inputs is inappropriate! (only even numbers)')
end


%for all input combinations
for iInput=1:totalNrOfInputs
    
    %initialize waitbar  
    %h = waitbar(0,['iInput:',num2str(iInput)]);    
    
    %generate input intervals
    uZ=segmentZonotope(inputField,iInput);
    
    %for all velocities
    for iVel=1:totalNrOfVelocities
            
        %obain current discrete state 
        iState=(iVel-1)*totalNrOfPositions+1;
        %update discretized state space
        stateField=set(stateField,'actualSegmentNr',iState);
        MC=set(MC,'field',stateField);
        
        %display iInput and iState
        iInput
        iState     
        
        %simulate hybrid automaton
        options.R0=segmentZonotope(stateField,iState);
        
        initialStates=gridPoints(interval(options.R0),5);
        sampleInputs=gridPoints(interval(uZ),5);
        
        %init finalStateMat
        finalStateMat.T=[];
        finalStateMat.OT=[];
        
        for initIndx=1:length(initialStates(1,:))
            %set initial state for simulation
            options.x0=initialStates(:,initIndx);
            for inputIndx=1:length(sampleInputs(1,:))
                %set input value
                for i=1:5
                    options.uLocTrans{i}=sampleInputs(:,inputIndx);
                    options.Uloc{i}=zonotope(0);
                end
                
                %determine initial location
                if (sampleInputs(:,inputIndx)>0) && (options.x0(2)<changeSpeed)
                    options.startLoc=1; %acceleration at slow speed
                elseif (sampleInputs(:,inputIndx)>0) && (options.x0(2)>=changeSpeed)
                    options.startLoc=2; %acceleration at high speed                 
                else
                    options.startLoc=4; %deceleartion
                end                  
                
                %set time steps
                for iTime=1:timeSteps
                    %set final time
                    options.tFinal=iTime*finalTime/timeSteps;
                    %simulate HA
                    HA=simulate(HA,options); 
                    %get final state
                    finalState=get(HA,'finalState');
                    %store result for time intervals
                    finalStateMat.OT(end+1,:)=finalState.x;
                end
                %store result for time points
                finalStateMat.T(end+1,:)=finalState.x;
                %get trajectory
                trajectories{initIndx,inputIndx}=get(HA,'trajectory');
            end
        end
        
        %Update Markov Chain
        MC=build4road(MC,finalStateMat,iInput,iState);
          %if mod(iState,19)==0
           if iState>0
                %plot(MC,trajectories,options,iState);
           end
        
        %update waitbar
        %waitbar(iState/totalNrOfStates,h); 
    end
end
%close waitbar
%close(h);

%optimize structure of Markov chain
[T, projMat, GammaFull] = convertTransitionMatrix(MC, gamma);

%create structure for abstract model
probModel.T = T;
probModel.projMat = projMat;
probModel.GammaFull = GammaFull;
probModel.timeStep = finalTime;
probModel.stateField = stateField;
probModel.inputField = inputField;


%------------- END OF CODE --------------
