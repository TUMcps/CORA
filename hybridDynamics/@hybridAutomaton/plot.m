function p=plot(obj,type,options)
% plot - plots results of a hybrid automaton;
% currently: simulation or reachability results
%
% Syntax:  
%    plot(obj,type)
%
% Inputs:
%    obj - hybrid automaton object
%    type - what should be plotted? 'simulation' or 'reachableSet'
%
% Outputs:
%    p - probability of hitting unsafe set for probabilistic plotting
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      11-May-2007 
% Last update:  28-September-2007
%               03-September-2009
% Last revision:---

%------------- BEGIN CODE --------------

%choice of plot type
p=[];
switch type
    case 'simulation'
        plotSim(obj,options);
    case 'simulationFromOrigin'
        plotSimFromOrigin(obj,options);
    case 'reachableSet'
        plotReach(obj,options);
    case 'reachableSetFilled'
        plotReachFilled(obj,options);
    case 'yawSet'
        plotYaw(obj,options);
    case 'reachableSetIntersections'
        plotReachIntersections(obj,options);        
    case 'probReachableSet'
        p=plotReachProb(obj,options);    
    case 'carInteraction'
        plotCarInteraction(obj);          
    otherwise
        disp('wrong plot type');
end



%------------- END OF CODE --------------