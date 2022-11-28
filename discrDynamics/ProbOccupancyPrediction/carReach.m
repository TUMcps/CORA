function probModel = carReach(fileName,pathName,modelInitialization)
% carReach - generates Markov model based on simulation
%
% Syntax:  
%    carReach(FileName,PathName)
%
% Inputs:
%    fileName - file name of fArray
%    pathName - path name of fArray
%    modelInitialization - handle to model initialization
%
% Outputs:
%    probModel - probabilistic model of a vehicle following a road
%
% Example: 
%
% 
% Author:       Matthias Althoff
% Written:      31-July-2016
% Last update:  31-July-2017
% Last revision:---


%------------- BEGIN CODE --------------

%set path
global filePath
filePath = [CORAROOT filesep 'contDynamics' filesep 'stateSpaceModels'];

%load fArray to determine segment length of road 
cd(pathName);
file=load(fileName);
fArray=file.fArray;

%load car model
[HA,options,params,stateField,inputField,changeSpeed] = modelInitialization(fArray.segLengthOther);

%set gamma value (how often inputs change)
gamma = 0.2;

%set final time and time steps
finalTime=params.tFinal;
timeSteps=5;

%Initialize Markov Chain
MC=markovchain(stateField);

%obtain number of segments of the state discretization
nrOfSegments = stateField.nrOfSegments;

%total number of discrete inputs, states, positions and velocities
totalNrOfInputs = nrOfCells(inputField);
totalNrOfPositions=nrOfSegments(1);
totalNrOfVelocities=nrOfSegments(2);

% check for even number of inputs
if rem(totalNrOfInputs,2)>0
        disp('Number of Inputs is inappropriate! (only even numbers)')
end


%for all input combinations
for iInput=1:totalNrOfInputs   
    
    %generate input intervals
    aux=cellZonotopes(inputField,iInput);
    uZ=aux{1};
    
    %for all velocities
    for iVel=1:totalNrOfVelocities
            
        %obain current discrete state 
        iState=(iVel-1)*totalNrOfPositions+1;
        
        %display iInput and iState
        iInput
        iState     
        
        %obtain intial set, sate samples, and inputv samples
        aux=cellZonotopes(stateField,iState);
        params.R0=aux{1};
        initialStates=gridPoints(interval(params.R0),5);
        sampleInputs=gridPoints(interval(uZ),5);
        
        %init finalStateMat
        finalStateMat.T=[];
        finalStateMat.OT=[];
        
        for initIndx=1:length(initialStates(1,:))
            %set initial state for simulation
            params.x0=initialStates(:,initIndx);
            for inputIndx=1:length(sampleInputs(1,:))
                %set input value
                for i=1:timeSteps
                    params.uLocTrans{i}=sampleInputs(:,inputIndx);
                    params.Uloc{i}=zonotope(0);
                end
                
                %determine initial location
                if (sampleInputs(:,inputIndx)>0) && (params.x0(2)<changeSpeed)
                    params.startLoc=1; %acceleration at slow speed
                elseif (sampleInputs(:,inputIndx)>0) && (params.x0(2)>=changeSpeed)
                    params.startLoc=2; %acceleration at high speed                 
                else
                    params.startLoc=4; %deceleartion
                end                  
                
                traj.t = {}; traj.x = {}; traj.loc = {};
                
                %set time steps
                for iTime=1:timeSteps
                    %set final time
                    params.tFinal=iTime*finalTime/timeSteps;
                    %simulate HA
                    [t,x,loc] = simulate(HA,params); 
                    traj.t = [traj.t;t];
                    traj.x = [traj.x;x];
                    traj.loc = [traj.loc;loc];
                    %get final state
                    finalState = x{end}(end,:)';
                    %store result for time intervals
                    finalStateMat.OT(end+1,:)=finalState;
                end
                %store result for time points
                finalStateMat.T(end+1,:)=finalState;
                %get trajectory
                trajectories{initIndx,inputIndx} = traj;
            end
        end
        
        %Update Markov Chain
        MC=build4road(MC,finalStateMat,iInput,iState);
          %if mod(iState,19)==0
%            if iState>0
%                 plot(MC,trajectories,options,iState);
%            end
    end
end

%optimize structure of Markov chain
[T, projMat, GammaFull] = convertTransitionMatrix(MC, gamma);

%create structure for abstract model
probModel.T = T;
probModel.projMat = projMat;
probModel.GammaFull = GammaFull;
probModel.timeStep = finalTime;
probModel.stateField = stateField;
probModel.inputField = inputField;
probModel.fArray = fArray;



%------------- END OF CODE --------------