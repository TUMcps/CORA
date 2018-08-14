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

%set path
global filePath
filePath = [coraroot '/contDynamics/stateSpaceModels'];

%load fArray to determine segment length of road 
cd(pathName);
file=load(fileName);
fArray=file.fArray;

%load car model
[HA,options,stateField,inputField] = modelInitialization(fArray.segLengthOther);

%set gamma value (how often inputs change)
gamma = 0.2;

%set final time and time steps
finalTime=options.tFinal;
 
%Initialize Markov Chain
MC=markovchain(stateField);

%obtain number of segments of the state discretization
%nrOfSegments=get(stateField,'nrOfSegments'); %<--AP
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
    
    %initialize waitbar  
    %h = waitbar(0,['iInput:',num2str(iInput)]);    
    
    %generate input intervals
    uZCell=cellZonotopes(inputField,iInput);
    uZ = uZCell{1}
    for i=1:length(get(HA,'location'))
        options.uLocTrans{i} = center(uZ);
        options.Uloc{i} = uZ + (-options.uLocTrans{i});
        options.uLoc{i} = center(uZ);
        options.timeStepLoc{i} = options.timeStep;
    end
    
    
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
        options.R0 = R0Cell{1};
        
        %determine initial location
        if iInput>(totalNrOfInputs/2)
            options.startLoc=1; %acceleration
        else
            options.startLoc=3; %deceleartion
        end           
        
        %compute reachable set
        HA=reach(HA,options);     
   
        %Update Markov Chain
        MC=build4road_reach(MC,HA,iInput,iState);
          %if mod(iState,19)==0
          %if iState>3000
                %plot_reach(MC,HA,options,iState);
          %end
        
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
probModel.fArray = fArray;

%------------- END OF CODE --------------
