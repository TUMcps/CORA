function [inputProb] = speedRestriction(simOptions,markovChainSpec)
% speedRestriction - Restricts the speed on roads due to speed limits or
% the road geometry (maximum speed through curves)
%
% Syntax:  
%    [inputProb] = speedRestriction(simOptions,markovChainSpec)
%    
%
% Inputs:
%    simOptions - simulation options
%    markovChainSpec - Markov-Chain specification
%
% Outputs:
%    inputProb - input probability vector
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      03-July-2008 
% Last update:  28-August-2008
% Last revision: ---

%------------- BEGIN CODE --------------

%obtain variables
profileHandle=simOptions.profileHandle;
stateField=simOptions.stateField;
inputField=simOptions.inputField;
%qFree=simOptions.freeDrivingProb;
T=markovChainSpec.timeStep;

%obtain acceleration values to parameterize the speed profiles
accVector(1)=input2acceleration_old(-1,NaN,simOptions.type); %lowest acceleration
iInput=2;
while accVector(iInput-1)<0
    %obtain left limit
    interval=cellIntervals(inputField,iInput);
    accVector(iInput)=input2acceleration_old(interval(1),NaN,simOptions.type); %lowest acceleration
    iInput=iInput+1;    
end
accVector(end)=[];

%get number of state and input segments
nrOfInputs=prod(get(inputField,'nrOfSegments'));
nrOfStates=prod(get(stateField,'nrOfSegments'));

%loop for maximum acceleartion
for iMaxAcc=1:length(accVector)
    %consider outside cell
    nrOfCrashs=nrOfInputs-1;
    inputProbTmp(iMaxAcc,1,:)=[1,zeros(1,nrOfCrashs)]; 
    
    %loop for each cell
    for iState=1:nrOfStates
        %get minimum position and velocity of the current cell
        interval=cellIntervals(stateField,iState);
        pos=interval(1,1);
        vel=interval(2,1);
        
        %get velocity of the velocity profile
        velProfile=profileHandle(pos,accVector(iMaxAcc));
        
        %init number of crashes (crash with vehicle that moves along
        %the velocity profile)
        nrOfCrashs=0;
        %loop for different accelerations
        for iInput=1:nrOfInputs
            %obtain acceleration
            interval=cellIntervals(inputField,iInput);
            acc=input2acceleration_old(interval(1),vel,simOptions.type); %lowest acceleration

            %check if the velocity is below the profile after time T
            newVel=acc*T+vel;
            if newVel>0
                newPos=0.5*acc*T^2+vel*T+pos;

                newVelProfile=profileHandle(newPos,accVector(iMaxAcc));

                %check if crash has happened
                if newVel>newVelProfile
                    nrOfCrashs=nrOfCrashs+1;
                end
            end
        end
        
        %maximum input for the vehicle
        safeInputs=nrOfInputs-nrOfCrashs;
        if safeInputs==0
            safeInputs=1;
            nrOfCrashs=nrOfInputs-safeInputs;
        end             
        
        %generate restriction vector
        resVec=[ones(1,safeInputs),zeros(1,nrOfCrashs)]; %no alpha values 
        inputProbTmp(iMaxAcc,iState+1,:)=resVec;      
        
    end
end

%combine input probabilities
accProb=1/length(accVector)*ones(length(accVector),1);
inputProb=0*squeeze(inputProbTmp(1,:,:));
for iMaxAcc=1:length(accVector)
    inputProb=inputProb+accProb(iMaxAcc)*squeeze(inputProbTmp(iMaxAcc,:,:));
end


%------------- END OF CODE --------------