function res = test_intersectionDatabase(~)
% test_intersectionDatabase - unit test function for testing the correct
% computation of the intersection database for probabilistic occupancy
% prediction of road vehicles
%
% Syntax:  
%    res = test_intersectionDatabase(~)
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      01-August-2017 
% Last update:  
% Last revision: ---

%------------- BEGIN CODE --------------


% obtain uncertain center intervalhull
% other traffic participants
parameters.segLength=4; %[m] 
parameters.centerWidth=3.5; %[m]
parameters.devSegments=3; 

% ego vehicle
parameters.segLengthEgo=1; %[m] 
parameters.segWidthEgo=0.5; %[m]

% obtain uncertain body interval
% vehicle data
parameters.carLength=4; %[m] 
parameters.carWidth=2; %[m]

parameters.bicycleLength=2; %[m] 
parameters.bicycleWidth=0.5; %[m]

% set number of x, y and angle segments
parameters.nrOfxSeg = 3;
parameters.nrOfySeg = 3;
parameters.nrOfAngleSeg = 3;

% compute intersection database
fArray = intersectionDatabase(parameters);

% load ground-truth results------------------------------------------------
load fArray_groundTruth
%--------------------------------------------------------------------------

% perform comparisons------------------------------------------------------
res_partial = zeros(length(fArray.val.car) + length(fArray.val.bicycle),1);
for iAngleSeg = 1:length(fArray.val.car)
    res_partial(iAngleSeg) = (max(max(abs(fArray.val.car{iAngleSeg} - fArray_groundTruth.val.car{iAngleSeg}))) < 1e-12);
end
for iAngleSeg = 1:length(fArray.val.bicycle)
    res_partial(iAngleSeg + length(fArray.val.car)) = (max(max(abs(fArray.val.bicycle{iAngleSeg} - fArray_groundTruth.val.bicycle{iAngleSeg}))) < 1e-12);
end

% have all partial tests passed?
res = prod(res_partial);
%--------------------------------------------------------------------------



%------------- END OF CODE --------------