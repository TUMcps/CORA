function res = testSpecial_markovchain_carReach
% testSpecial_markovchain_carReach - unit_test_function for building a
% Markov chain of a car using simulation results
%
% Syntax:
%    res = testSpecial_markovchain_carReach
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false

% Authors:       Matthias Althoff
% Written:       31-July-2016
% Last update:   07-July-2021
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%set path
global filePath
filePath = [CORAROOT filesep 'contDynamics' filesep 'stateSpaceModels'];

%load car model
[HA,params,params,stateField,inputField,changeSpeed] = testAuxiliaryFct_initCar;

%set gamma value (how often inputs change)
gamma = 0.2;

%set final time and time steps
finalTime = 0.5;
timeSteps = 5;

%Initialize Markov Chain
MC = markovchain(stateField);

%obtain number of segments of the state discretization
nrOfSegments = stateField.nrOfSegments;

%total number of discrete inputs, states, positions and velocities
totalNrOfInputs = prod(inputField.nrOfSegments);
totalNrOfPositions = nrOfSegments(1);
totalNrOfVelocities = nrOfSegments(2);

% check for even number of inputs
if rem(totalNrOfInputs,2)>0
        disp('Number of Inputs is inappropriate! (only even numbers)')
end


%for all input combinations
for iInput = 1:totalNrOfInputs 
    
    %generate input intervals
    tmp = cellZonotopes(inputField,iInput);
    uZ = tmp{1};
    
    %for all velocities
    for iVel = 1:totalNrOfVelocities
            
        %obain current discrete state 
        iState = (iVel-1)*totalNrOfPositions+1;
%         %update discretized state space
%         stateField = set(stateField,'actualSegmentNr',iState);
%         MC = set(MC,'field',stateField);
        
        %display iInput and iState
        iInput
        iState     
        
        %simulate hybrid automaton
        tmp = cellZonotopes(stateField,iState);
        params.R0 = tmp{1};
        
        initialStates = gridPoints(interval(params.R0),5);
        sampleInputs = gridPoints(interval(uZ),5);
        
        %init finalStateMat
        finalStateMat.T = [];
        finalStateMat.OT = [];
        
        for initIndx = 1:length(initialStates(1,:))
            %set initial state for simulation
            params.x0 = initialStates(:,initIndx);
            for inputIndx = 1:length(sampleInputs(1,:))
                %set input value
                for i = 1:5
                    params.uLocTrans{i} = sampleInputs(:,inputIndx);
                    params.Uloc{i} = zonotope(0);
                end
                
                %determine initial location
                if (sampleInputs(:,inputIndx)>0) && (params.x0(2)<changeSpeed)
                    params.startLoc = 1; %acceleration at slow speed
                elseif (sampleInputs(:,inputIndx)>0) && (params.x0(2)>=changeSpeed)
                    params.startLoc = 2; %acceleration at high speed                 
                else
                    params.startLoc = 4; %deceleartion
                end                  
                
                %set time steps
                for iTime = 1:timeSteps
                    %set final time
                    params.tFinal = iTime*finalTime/timeSteps;
                    %simulate HA
                    [~,x] = simulate(HA,params); 
                    %get final state
                    finalState = x{end}(end,:);
                    %store result for time intervals
                    finalStateMat.OT(end+1,:) = finalState;
                end
                %store result for time points
                finalStateMat.T(end+1,:) = finalState;
                %get trajectory
                trajectories{initIndx,inputIndx} = x;
            end
        end
        
        %Update Markov Chain
        MC = build4road(MC,finalStateMat,iInput,iState);
          %if mod(iState,19)==0
           if iState>0
                %plot(MC,trajectories,params,iState);
           end
        
    end
end

%optimize structure of Markov chain
[T, projMat, GammaFull] = convertTransitionMatrix(MC, gamma);

%% compute maximum errors
% init
resPartial = [];
accuracy = 1e-8;

% load precomputed results
load carReach_unitTest Tsave projMatSave GammaFullSave

% transition matrix - discrete time
error = abs(T.T - Tsave.T);
maxError = max(max(error));
resPartial(end+1) = (maxError < accuracy);

% transition matrix - continuous time
error = abs(T.OT - Tsave.OT);
maxError = max(max(error));
resPartial(end+1) = (maxError < accuracy);

% projection matrix
error = abs(projMat - projMatSave);
maxError = max(max(error));
resPartial(end+1) = (maxError < accuracy);

% Gamma matrix
error = abs(GammaFull - GammaFullSave);
maxError = max(max(error));
resPartial(end+1) = (maxError < accuracy);

% final result
res = all(resPartial);

% ------------------------------ END OF CODE ------------------------------
