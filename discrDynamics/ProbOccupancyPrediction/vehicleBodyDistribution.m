function [probVec]=vehicleBodyDistribution(oldPartition,newPartition,stretch,probVec)
% vehicleBodyDistribution - transforms the center point probability
% distribution of a vehicle to the distribution of the vehicle body in
% lateral or longitudinal direction; this can also be done for a new
% partition of the obtained distribution
%
% Syntax:  
%    [probVec]=vehicleBodyDistribution(oldPartition,newPartition,stretch,pr
%    obVec)
%
% Inputs:
%    oldPartition - old partition object
%    newPartition - new partition object
%    stretch - lateral or longitudinal length of the vehicle
%    probVec - probability vector of the partition
%
% Outputs:
%    probVec - new probability vector of the partition
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      27-January-2009
% Last update:  25-July-2016
%               13-November-2017
% Last revision:---

%------------- BEGIN CODE ---------------

%determine non-zero cells
cellInd=find(probVec);

%set up interval for the stretching
IHadd=interval(-0.5*stretch,0.5*stretch);

%initialize total probability vector
nrOfSegments=newPartition.nrOfSegments; %number of segments of new partition
tpTotal=zeros(nrOfSegments+1,1); %0 for the outside probability

%compute projections of each partition interval onto the new partition
for i=1:length(cellInd)
    %obtain interval of selected cell
    int=cellIntervals(oldPartition,cellInd(i));
    IH=int{1};
    
    %stertch interval
    IH=IH+IHadd;
    
    %intersect interval with new partition
    [~, iProp] = exactIntersectingCells(newPartition,IH);
    tpTmp=probVec(cellInd(i))*iProp;
    
    %sum up to total probability
    tpTotal=tpTotal+tpTmp;
end

%overwrite probability vector
probVec=tpTotal(2:end); %outside probability "thrown away"

%------------- END OF CODE --------------