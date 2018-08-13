function [options] = timeStep(options,Rnext)
% timeStep - adjusts the time step of the reachable set computation
%
% Syntax:  
%    [options] = timeStep(options)
%
% Inputs:
%    options - options struct
%
% Outputs:
%    options - options struct
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: 

% Author: Matthias Althoff
% Written: 15-January-2008
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------

explorationSpeed=abs(Rnext.f0max);

%time step is the time that maximum speed needs to travel weighted with
%explorationWeight
options.timeStep=min(options.explorationWeight./explorationSpeed');

%check if time boundary is respected
if options.timeStep<options.timeStepInterval(1)
    options.timeStep=options.timeStepInterval(1);
elseif options.timeStep>options.timeStepInterval(2)
    options.timeStep=options.timeStepInterval(2);
end

%------------- END OF CODE --------------