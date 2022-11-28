function probModel = carReach_reach(fileName,pathName,modelInitialization)
% carReach_reach - generates Markov model based on reachability analysis
%
% Syntax:  
%    probModel = carReach_reach(fileName,pathName,modelInitialization)
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

%load fArray to determine segment length of road 
cd(CORAROOT);
file=load(fileName);
fArray=file.fArray;

%load car model
[HA,options,params,stateField,inputField] = modelInitialization(fArray.segLengthOther);

%set gamma value (how often inputs change)
gamma = 0.2;

%set final time and time steps
finalTime=params.tFinal;
 
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
        disp('Number of Inputs is inappropriate! (only even number)')
end


%for all input combinations
for iInput=1:totalNrOfInputs  
    
    %generate input intervals
    uZCell=cellZonotopes(inputField,iInput);
    uZ = uZCell{1};
    params.U = uZ;
    
    %for all velocities
    for iVel=1:totalNrOfVelocities
            
        %obain current discrete state 
        iState=(iVel-1)*totalNrOfPositions+1;
        %update discretized state space
        MC=set(MC,'field',stateField);
        
        %display iInput and iState
        iInput
        iState
        
        R0Cell=cellZonotopes(stateField,iState);
        params.R0 = R0Cell{1};
        
        %determine initial location
        if iInput>(totalNrOfInputs/2)
            params.startLoc=1; %acceleration
        else
            params.startLoc=3; %deceleartion
        end           
        
        %compute reachable set
        R = reach(HA,params,options);     
   
        %Update Markov Chain
        MC=build4road_reach(MC,R,iInput,iState);
%           if mod(iState,19)==0
%           if iState>3000
%                 plot_reach(MC,HA,R,params,iState);
%           end
        
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
