function [pChange]=potChange(simOptions)
% potChange - determines the conditional probability that starting from
% a certain discrete state, one can change the lane without violating the 
% vehicle that will be the following vehicle after the lane change
%
% Syntax:  
%    [pChange]=potChange(simOptions)
%
% Inputs:
%    simOptions - simulation options
%
% Outputs:
%    pChange - probability vector of the discrete states that the lane
%    change can be performed due to the following vehicle
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author: Matthias Althoff
% Written: 19-July-2008 
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------

%obtain number of segments
field=simOptions.field;
nrOfSeg=get(field,'nrOfSegments');
nrOfCells=prod(nrOfSeg);

%set iMode for which the following car is not forced to break
iMode=floor(simOptions.nrOfInputs/2)+1; %use probability distribution for this?

for iStep=1:simOptions.runs
    for iCell=1:nrOfCells
        for iAheadMode=1:simOptions.nrOfInputs
            %obtain 1/0 selection matrix mBool
            mBool=find(Omega{iMode,iAheadMode}(:,iCell));
            %obtain change probability
            pChange{iStep}(iCell)=sum(mBool.*simOptions.leftBehindProbVector{iStep});
        end
    end
end


%------------- END OF CODE --------------