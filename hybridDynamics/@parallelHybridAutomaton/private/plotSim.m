function plotSim(obj, options)
% plotSim - plot the simulation results stored in a parallel hybrid 
%           automaton object 
%
% Syntax:  
%    plotSim(obj,options)
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
% Written:      04-July-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

dim = options.projectedDimensions;
hold on

% loop over all stored simulations
for i = 1:length(obj.result.simulation)
   
    simRes = obj.result.simulation{i};
    
    % loop over all visited location
    for j = 1:length(simRes.x)
        
       plot(simRes.x{j}(:,dim(1)),simRes.x{j}(:,dim(2)),options.plotType);
    end  
end

%------------- END OF CODE --------------