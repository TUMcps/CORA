function R = reach(obj,options)
% reach - computes the reachable sets of the discrete time system
%
% Syntax:  
%    R = reach(obj,options)
%
% Inputs:
%    obj - nonlinearSysDT object
%    options - options for the computation of the reachable set
%
% Outputs:
%    R - cell array containing the reachable sets
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff, Niklas Kochdumper
% Written:      21-August-2012
% Last update:  29-January-2018
% Last revision:---

%------------- BEGIN CODE --------------


% initialize cell array that stores the reachable sets
steps = (options.tFinal-options.tStart)/options.timeStep;
R = cell(steps+1,1);
R{1} = options.R0;

t=options.tStart;
iSet=1;

while t<options.tFinal
    
    % if a trajectory should be tracked
    if isfield(options,'uTransVec')
        options.uTrans = options.uTransVec(:,iSet);
    end  
    
    % compute next reachable set
    R{iSet+1} = linReach(obj,R{iSet},options);
    
    % increment time and set counter
    t = t+options.timeStep;
    iSet = iSet+1; 
    options.t=t;
    if isfield(options,'verbose') && options.verbose 
        disp(t); %plot time
    end
end


%------------- END OF CODE --------------