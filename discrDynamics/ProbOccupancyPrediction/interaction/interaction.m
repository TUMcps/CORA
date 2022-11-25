function ThetaC = interaction(fileName,pathName,modelInitialization)
% interaction - computes interaction matrix Theta for two road vehicles
% specified by the same dynamic model
%
% Syntax:  
%    
%    ThetaC = interaction(fileName,pathName,modelInitialization)
%
% Inputs:
%    fileName - file name of fArray
%    pathName - path name of fArray
%    modelInitialization - handle to model initialization
%
% Outputs:
%    ThetaC - interaction matrix (optimized structure)
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
% Written:      28-March-2008 
% Last update:  27-June-2008
%               10-August-2018
%               12-October-2009
%               10-August-2018
% Last revision: ---

%------------- BEGIN CODE --------------

%set path
global filePath
filePath = [coraroot '/contDynamics/stateSpaceModels'];

%load fArray to determine segment length of road 
cd(pathName);
file=load(fileName);
fArray=file.fArray;

%load car model
[~,options,stateField,inputField] = modelInitialization(fArray.segLengthOther);

%define attention value epsilon
epsilon=1e-4;

%define car length
carLength=5;

%set time increment
T=options.tFinal;

%car F (following) and car L (leading) have the same state partitions than 
%specified for the Markov chain generation
carFstates=stateField;
carLstates=stateField;

%car F and car L have the same input partitions than specified for the
%Markov chain generation
carFinputs=inputField;
carLinputs=inputField;

%number of states
nrOfCarFstates=nrOfCells(carFstates);
nrOfCarLstates=nrOfCells(carLstates);

%number of inputs
nrOfInputs=nrOfCells(inputField);

%initialize entries
for carLmode=1:nrOfInputs
    for carFmode=1:nrOfInputs
        Theta{carFmode,carLmode}=sparse(nrOfCarFstates+1,nrOfCarLstates+1);
    end
end

%check possible state interval combinations
%car F is behind car L
for carFcell=1:nrOfCarFstates
    for carLcell=1:nrOfCarLstates   
        %track progress
        carFcell
        carLcell
        
        %get cell intervals of car F
        c = cellCenter(carFstates,carFcell);
        carFstate = c{1};
        carFstate(1)=carFstate(1)+carLength; %consider car length
        %get cell intervals of car L
        c = cellCenter(carLstates,carLcell);
        carLstate = c{1};
        %car F behind car L?
        FbehindL=(carFstate(1)<carLstate(1));
        
        if FbehindL
            %combination of different acceleration intervals
            for carLmode=1:nrOfInputs
                
                %init deltaTVec
                deltaTVec=[T 4*T 8*T];
                
                for iDeltaT=1:3
                
                    %obtain deltaT
                    deltaT=deltaTVec(iDeltaT);

                    %init carFmode, crash
                    carFmode=0;
                    crash=0;
                    safeInputs(iDeltaT)=0;

                    while (crash==0) && (carFmode<nrOfInputs)
                        %increase carFmode
                        carFmode=carFmode+1;
      
                        %get input intervals of cars
                        c = cellCenter(carFinputs,carFmode);   
                        carFinput = c{1};
                        c = cellCenter(carLinputs,carLmode);
                        carLinput = c{1};

                        %compute if a crash occurs
                        crash=crashCheck(carFstate,carLstate,carFinput,carLinput,deltaT);
                        safeInputs(iDeltaT)=safeInputs(iDeltaT)+(1-crash);
                    end
                    %obtain number of crashes
                    nrOfCrashs(iDeltaT)=nrOfInputs-safeInputs(iDeltaT);

                    if safeInputs(iDeltaT)==0
                        safeInputs(iDeltaT)=1;
                        nrOfCrashs(iDeltaT)=nrOfInputs-safeInputs(iDeltaT);
                    end

                    %generate restriction vector
                    resVec(iDeltaT,:)=[ones(1,safeInputs(iDeltaT)),epsilon*ones(1,nrOfCrashs(iDeltaT))];
                end
                
               
                %compute overall restriction vector
                resVecAvg=1/3*(resVec(1,:)+resVec(2,:)+resVec(3,:));

                %update Theta
                for carFmode=1:nrOfInputs
                    Theta{carFmode,carLmode}(carFcell+1,carLcell+1)=sparse(resVecAvg(carFmode));
                end
            end  
        else
            %car F is in front of car L
            for carLmode=1:nrOfInputs
                safeInputs=1;
                nrOfCrashs=nrOfInputs-safeInputs;  
                resVecAvg=[ones(1,safeInputs),epsilon*ones(1,nrOfCrashs)];
                %update Theta
                for carFmode=1:nrOfInputs
                    Theta{carFmode,carLmode}(carFcell+1,carLcell+1)=sparse(resVecAvg(carFmode));
                end      
            end
        end
    end  
end

%set Theta values for the outside state
for carLmode=1:nrOfInputs
    Theta{1,carLmode}(1,:)=ones(1,nrOfCarFstates+1);
end

% obtain optimized structure of interaction matrix
[ThetaC]=convertInteractionMatrix(Theta);


%------------- END OF CODE --------------