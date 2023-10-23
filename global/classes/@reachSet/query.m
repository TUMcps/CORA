function val = query(R,prop)
% query - get properties of the reachable set
%
% Syntax:
%    val = query(R,prop)
%
% Inputs:
%    R - reachSet object
%    prop - property: 
%           'reachSet': cell-array of time-interval solutions
%           'reachSetTimePoint': cell-array of time-point solutions
%           'finalSet': final time-point solution
%           'tVec': vector of time steps
%           'tFinal': final time stamp
%           'allLoc': all location IDs of visited locations
%
% Outputs:
%    val - value of the property
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: reachSet

% Authors:       Niklas Kochdumper
% Written:       02-June-2020
% Last update:   19-May-2023 (MW, add 'allLoc' property)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % check input arguments
    inputArgsCheck({{R,'att','reachSet'}; ...
                    {prop,'str',{'reachSet','reachSetTimePoint',...
                    'finalSet','tVec','tFinal','allLoc'}}});

    switch prop
        case 'reachSet'
    
            val = R(1,1).timeInterval.set;
            for i = 2:size(R,1)
               val = [val;R(i).timeInterval.set]; 
            end
        
        case 'reachSetTimePoint'

            val = R(1,1).timePoint.set;
            for i = 2:size(R,1)
               val = [val;R(i).timePoint.set]; 
            end
        
        case 'finalSet'

            if size(R,1) == 1
               val = R(1,1).timePoint.set{end}; 
            else
               parents = aux_getParents(R);
               ind = setdiff(1:size(R,1),parents);
               val = cell(length(ind),1);
               
               for i = 1:length(ind)
                  val{i} = R(ind(i),1).timePoint.set{end};
               end
            end
        
        case 'tVec'
    
            if size(R,1) == 1
                val = diff(cell2mat(R.timePoint.time));
            else
                throw(CORAerror('CORA:notSupported',...
                    'Multiple branches not supported.'));
            end

        case 'tFinal'

            val = R(1,1).timePoint.time{end};
            for i = 2:size(R,1)
                val = max(val,R(i).timePoint.time{end});
            end
    
            if isa(val, 'interval')
                val = supremum(val);
            end

        case 'allLoc'

            % at least one branch in reachSet object
            val = R(1).loc;
            
            % loop over all reachable sets, find all location IDs
            for i=1:length(R)
                % check if location ID has been seen before
                if ~compareMatrices(R(i).loc,val,eps,'subset')
                    val = [val,R(i).loc];
                end
            end
         
    end
end


% Auxiliary functions -----------------------------------------------------

function parents = aux_getParents(obj)
% get the indices of all parent reachable sets

    parents = [];
    
    for i = 1:size(obj,1)
        parents = [parents;obj(i,1).parent];
    end

    parents = unique(parents);
end

% ------------------------------ END OF CODE ------------------------------
