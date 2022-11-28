function start_intersectionDatabase(~)
% start_intersectionDatabase - sets parameters and stores the results of
% the intersectionDatabase
%
% Syntax:  
%    start_intersectionDatabase(~)
%
% Inputs:
%    -
%
% Outputs:
%    -
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


% obtain uncertain center interval
% other traffic participants
parameters.segLength=4; %[m] 
parameters.centerWidth=3.5; %[m]
parameters.devSegments=7; 

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
parameters.nrOfxSeg = 44;
parameters.nrOfySeg = 44;
parameters.nrOfAngleSeg = 72;

% compute intersection database
fArray = intersectionDatabase(parameters);


% save intersection database
[file,path] = uiputfile('*.mat','Save Intersection Database As');
cd(path);
save(file,'fArray');




%------------- END OF CODE --------------