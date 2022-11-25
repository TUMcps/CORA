function [devProbLeft,dispDevProbLeft,devProbRight,dispDevProbRight]=deviationProbability(obj,lcEvolProb)
% deviationProbability - computes the lateral deviation probability based
% on the lane change probability
%
% Syntax:  
%    [p]=deviationProbability(R,changeRatio)
%
% Inputs:
%    obj - road object 
%    lcEvolProb - probability distribution collection for changing lanes 
%
% Outputs:
%    devProb - deviation probability distribution
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author: Matthias Althoff
% Written: 26-March-2009
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------

%get number of deviation segments
nrOfDevSegments=obj.nrOfDevSegments;

%set deviation probability distribution
if nrOfDevSegments==8
    devProbCarStd=[0,0,1,2,2,1,0,0];
else
    disp('Further deviation probabilities need to be defined.')
end

%normalize deviation probabilities
devProbCarStd=devProbCarStd/sum(devProbCarStd);

%check if lcEvolProb is empty
if isscalar(lcEvolProb) %only the number of runs is specified
    
    %set up difference between deviation and deviation display
    deviationFieldCenter=partition([-2, 2],8);
    deviationFieldBody=partition([-2, 2],7);
    stretch=2; %[m] 
    
    for i=1:lcEvolProb
        devProbLeft(i,:)=devProbCarStd;
        [dispDevProbLeft(i,:)]=vehicleBodyDistribution(deviationFieldCenter,deviationFieldBody,...
        stretch,devProbLeft(i,:));         
    end  
       
    devProbRight=[];
    dispDevProbRight=[];
else
    %compute deviations; 
    %assumptions: lane change takes 5 sec (10 steps), step size is 0.5 sec, 
    %road width is 4 m
    i=1;
    a=4/2.5^2;
    for j=1:5
        t=(j-1)*0.5;
        y(i)=0.5*a*t^2-2;
        i=i+1;
    end
    for j=1:6
        t=(j-1)*0.5;
        y(i)=4/2.5*t-0.5*a*t^2;
        i=i+1;
    end  
    
    %set up partition for the deviation intervals
    devField=partition([-4,4],...  %acceleartion in m/s^2
                     16);  
    devField2=partition([-2,2],...  %acceleartion in m/s^2
                     8);                  
    
    %generate interval hulls for deviation intervals
    [IH]=cellIntervals(devField2);
    

    for iStep=1:length(y)
        %translate interval hulls by the mean deviation y
        for j=1:length(IH)
            deltaIH{j}=IH{j}+y(iStep);
        end
        
        totalProb{iStep}=zeros(17,1);
        %intersect translated interval hulls with the cells and weight them
        %by the probabilities
        for j=1:length(IH)
            [~,prob]=exactIntersectingCells(devField,deltaIH{j});
            prob=prob*devProbCarStd(j);
            totalProb{iStep}=totalProb{iStep}+prob;
        end        
    end
    
    %set up difference between deviation and deviation display
    deviationFieldCenter=partition([-4, 4],16);
    deviationFieldBody=partition([-4, 4],14);
    stretch=2; %[m]    
       
    %compute probability distributions for each time step
    nrOfSteps=length(lcEvolProb(:,1));  
    
    for timeStep=1:nrOfSteps
        %init devProb
        devProb{timeStep}=0*totalProb{1};
        %sum results from different lc phases
        for iStep=1:length(totalProb)
            devProb{timeStep}=devProb{timeStep}+...
                lcEvolProb(timeStep,iStep)*totalProb{iStep};
        end
        
        %compute body distribution from vehicle center distribution
        [dispDevProb{timeStep}]=vehicleBodyDistribution(deviationFieldCenter,deviationFieldBody,...
        stretch,devProb{timeStep}(2:end));             
        
        %separate results for the left and right lane: center
        devProbLeft(timeStep,:)=devProb{timeStep}(2:9); 
        devProbRight(timeStep,:)=devProb{timeStep}(10:17); 
        
        %separate results for the left and right lane: body
        dispDevProbLeft(timeStep,:)=dispDevProb{timeStep}(1:7); 
        dispDevProbRight(timeStep,:)=dispDevProb{timeStep}(8:14);         
    end    
end

%------------- END OF CODE --------------