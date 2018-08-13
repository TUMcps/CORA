function plotYaw(obj,options)
% plot - plots reachable sets of a hybrid automaton
%
% Syntax:  
%    plotReach(obj)
%
% Inputs:
%    obj - hybrid automaton object
%
% Outputs:
%    none
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author: Matthias Althoff
% Written: 11-May-2007 
% Last update: 28-September-2007
%              26-March-2008
%              23-June-2009
% Last revision: ---

%------------- BEGIN CODE --------------


%load data from object structure
R=obj.result.reachSet.Rcont.OT;
loc=obj.result.reachSet.location;
dim=options.projectedDimensions;



%simplify sets
for (i=1:length(R))
    
    disp(['next set: ',num2str(i)]);
        
    %plot reachable sets of the same time interval
    for k=1:length(R{i}{1})
        %if reachable set  is a zonotope
        if strcmp('zonotope',class(R{i}{1}{k}))
            Rred=project(R{i}{1}{k},dim);
            Rred=reduce(Rred,'girard',options.polytopeOrder);
            plot(Rred,[1 2],options.plotType);
        else
            disp('set not plotted');
        end
    end
end




%------------- END OF CODE --------------