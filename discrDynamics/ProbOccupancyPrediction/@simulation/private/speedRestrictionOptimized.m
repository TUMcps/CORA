function [inputProb]=speedRestrictionOptimized(simOptions,markovChainSpec)
% speedRestriction - Restricts the speed on roads due to speed limits or
% the road geometry (maximum speed through curves)
%
% Syntax:  
%    [inputProb]=speedRestriction(simOptions,markovChainSpec)
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
%               01-July-2009
%               09-October-2009
%               13-November-2017
% Last revision: ---

%------------- BEGIN CODE --------------

%obtain variables
profileHandle=simOptions.profileHandle;
stateField=simOptions.stateField;
inputField=simOptions.inputField;
%qFree=simOptions.freeDrivingProb;
T=markovChainSpec.timeStep;

% %obtain acceleration values to parameterize the speed profiles
% accVector(1)=input2acceleration_old(-1,NaN,simOptions.type); %lowest acceleration
% iInput=2;
% while accVector(iInput-1)<0
%     %obtain left limit
%     interval=cellIntervals(inputField,iInput);
%     accVector(iInput)=input2acceleration_old(interval(1),NaN,simOptions.type); %lowest acceleration
%     iInput=iInput+1;    
% end
% accVector(end)=[];
accVector=linspace(7,1,4);

%get number of state and input segments
nrOfInputs=prod(inputField.nrOfSegments);
nrOfStates=prod(stateField.nrOfSegments);

%loop for maximum acceleartion
for iMaxAcc=1:length(accVector)
    %consider outside cell
    nrOfCrashs=nrOfInputs-1;
    inputProbTmp(1:nrOfInputs,iMaxAcc)=[1;zeros(nrOfCrashs,1)]; 
    
    %loop for each cell
    for iState=1:nrOfStates
        %get minimum position and velocity of the current cell
        interval=cellIntervals(stateField,iState);
        minValues = infimum(interval{1});
        pos=minValues(1);
        vel=minValues(2);
        
        %get velocity of the velocity profile
        velProfile=profileHandle(pos,accVector(iMaxAcc));
        
        %init number of crashes (crash with vehicle that moves along
        %the velocity profile)
        nrOfCrashs=0;
        %loop for different accelerations
        for iInput=1:nrOfInputs
            %obtain acceleration
            interval=cellIntervals(inputField,iInput);
            minValues = infimum(interval{1});
            acc=input2acceleration_old(minValues(1),vel,simOptions.type); %lowest acceleration

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
        resVec=[ones(safeInputs,1);zeros(nrOfCrashs,1)]; %no alpha values 
        newStates=nrOfInputs*iState+(1:nrOfInputs);
        inputProbTmp(newStates,iMaxAcc)=resVec;      
        
    end
end

%combine input probabilities
%accProb=1/length(accVector)*ones(length(accVector),1);
accProb=[0.05 0.05 0.1 0.8];
inputProb=0*inputProbTmp(:,1);
for iMaxAcc=1:length(accVector)
    inputProb=inputProb+accProb(iMaxAcc)*inputProbTmp(:,iMaxAcc);
end


%------------- END OF CODE --------------