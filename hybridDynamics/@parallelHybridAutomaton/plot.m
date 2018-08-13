function plot(obj,type,options)
% plot - plots results stored in a parallel hybrid automaton object
%        currently: simulation or reachability results
%
% Syntax:  
%    plot(obj,type,options)
%
% Inputs:
%    obj - hybrid automaton object
%    type - 'simulation' or 'reachableSet'
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
% Written:      04-July-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% choice of plot type
switch type
    case 'simulation'
        plotSim(obj,options);
    case 'reachableSet'
        plotReach(obj,options);      
    otherwise
        error('Wrong plot type!');
end

%------------- END OF CODE --------------