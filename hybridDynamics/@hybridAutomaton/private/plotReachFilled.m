function plotReachFilled(obj,options)
% plotReachFilled - plots reachable sets of a hybrid automaton
%
% Syntax:  
%    plotReachFilled(obj)
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

% Author: Niklas Kochdumper
% Written: 26-June-2018
% Last update: ---
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
                plotFilled(Rtemp,dimTemp,options.plotType{loc(i)},'EdgeColor','none');
            else
                plotFilled(Rtemp,dimTemp,options.plotType,'EdgeColor','none');
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
                    plotFilled(Rtemp,dimTemp,options.plotType{loc(i)},'EdgeColor','none');
                else
                    plotFilled(Rtemp,dimTemp,options.plotType,'EdgeColor','none');
                end
            end
        end
    end
end


%------------- END OF CODE --------------