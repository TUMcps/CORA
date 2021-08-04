function val = query(R,prop)
% query - get properties of the reachable set
%
% Syntax:  
%    val = query(R,prop)
%
% Inputs:
%    R - reachSet object
%    prop - property ('reachSet', 'reachSetTimePoint', or 'finalSet')
%
% Outputs:
%    val - value of the property
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: reachSet

% Author:       Niklas Kochdumper
% Written:      02-June-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    if strcmp(prop,'reachSet')
        val = R(1,1).timeInterval.set;
        for i = 2:size(R,1)
           val = [val;R(i).timeInterval.set]; 
        end
        
    elseif strcmp(prop,'reachSetTimePoint')
        val = R(1,1).timePoint.set;
        for i = 2:size(R,1)
           val = [val;R(i).timePoint.set]; 
        end
        
    elseif strcmp(prop,'finalSet')
        if size(R,1) == 1
           val = R(1,1).timePoint.set{end}; 
        else
           parents = getParents(R);
           ind = setdiff(1:size(R,1),parents);
           val = cell(length(ind),1);
           
           for i = 1:length(ind)
              val{i} = R(ind(i),1).timePoint.set{end};
           end
        end
        
    elseif strcmp(prop,'tVec')
        if size(R,1) == 1
            val = [R.timePoint.time{1}; diff(cell2mat(R.timePoint.time))];
        else
            error("Currently not implemented");
        end
        
    else
       error('Wrong value for input arguments "prop"!'); 
    end
end


% Auxiliary Functions -----------------------------------------------------

function parents = getParents(obj)
% get the indices of all parent reachable sets

    parents = [];
    
    for i = 1:size(obj,1)
        parents = [parents;obj(i,1).parent];
    end

    parents = unique(parents);
end

%------------- END OF CODE --------------