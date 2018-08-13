function plotReach(obj, options)
% plotReach - plot the reachable set stored in a parallel hybrid 
%             automaton object 
%
% Syntax:  
%    plotReach(obj,options)
%
% Inputs:
%    obj - parallel hybrid automaton object
%    options - structure with additional plot options
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

% Author:       Niklas Kochdumper
% Written:      05-July-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

dim = options.projectedDimensions;
hold on

% loop over all visited locations
for i = 1:length(obj.result.reachSet)
   
    reachRes = obj.result.reachSet{i};
    
    % loop over all reachable sets for the current location
    for j = 1:length(reachRes.R.OT)
        
       plot(reachRes.R.OT{j},dim,options.plotType);
    end  
end

%------------- END OF CODE --------------