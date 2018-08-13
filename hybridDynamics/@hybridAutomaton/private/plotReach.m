function plotReach(obj,options)
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
%              10-November-2010
%              23-November-2017
% Last revision: ---

%------------- BEGIN CODE --------------


%load data from object structure
R=obj.result.reachSet.R.OT; %<-- change back
%R=obj.result.reachSet.R.OT;
%R=obj.result.reachSet.R.OT{1};
loc=obj.result.reachSet.location;
dim=options.projectedDimensions;


for i=1:(length(R))
    disp(['next plot: ',num2str(i)]);
    for j=1:length(R{i})
        if ~iscell(R{i}{j}) % no split sets
            %plot reachable set
            Rtemp = R{i}{j};
            dimTemp = dim;
            if isa(Rtemp,'zonotope') || isa(Rtemp,'zonotopeBundle')
                Rtemp = project(Rtemp,dim);
                Rtemp = reduce(Rtemp,'girard',options.polytopeOrder);
                dimTemp = [1,2];
            end
            if iscell(options.plotType)
                plot(Rtemp,dimTemp,options.plotType{loc(i)});
            else
                plot(Rtemp,dimTemp,options.plotType);
            end
        else % sets are split
            for iSubSet=1:length(R{i}{j})
                %plot reachable set
                Rtemp = R{i}{j}{iSubSet};
                dimTemp = dim;
                if isa(Rtemp,'zonotope') || isa(Rtemp,'zonotopeBundle')
                    Rtemp = project(Rtemp,dim);
                    Rtemp = reduce(Rtemp,'girard',options.polytopeOrder);
                    dimTemp = [1,2];
                end
                if iscell(options.plotType)
                    plot(Rtemp,dimTemp,options.plotType{loc(i)});
                else
                    plot(Rtemp,dimTemp,options.plotType);
                end
            end
        end
    end
end


%------------- END OF CODE --------------