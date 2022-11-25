function res = find(R,prop,val)
% find - get reachSet object that satisfy a certain condition
%
% Syntax:  
%    res = find(R,prop,val)
%
% Inputs:
%    R - reachSet object
%    prop - property for condition ('location', 'parent', 'time')
%    val - value for property
%
% Outputs:
%    res - all reachSet objects that satisfy the condition
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

    res = [];
    
    % get all reachSet objects with the specified location
    if strcmp(prop,'location')
        
        for i = 1:size(R,1)
           if R(i).loc == val
              if isempty(res)
                  res = R(i); 
              else
                  res = add(res,R(i));
              end
           end
        end
        
    % get all reachSet objects with the specified parent
    elseif strcmp(prop,'parent')
        
        for i = 1:size(R,1)
           if R(i).parent == val
              if isempty(res)
                  res = R(i); 
              else
                  res = add(res,R(i));
              end
           end
        end
        
    % get all reachSet objects inside the specified time interval
    elseif strcmp(prop,'time')
        
        if ~isa(val,'interval')
            val = interval(val);
        end
        
        for i = 1:size(R,1)
            
           ind = [];
           
           for j = 1:length(R(i).timeInterval.time) 
              if isIntersecting(R(i).timeInterval.time{j},val)
                  ind = [ind;j];
              end
           end
            
           temp = getSubset(R(i),ind);
           
           if isempty(res)
               res = temp; 
           else
               res = add(res,temp);
           end   
        end
    end    
end


% Auxiliary Functions -----------------------------------------------------

function R = getSubset(R,ind)
% get the subset of the reachSet object corresponding to the indices

    timePoint.set = R.timePoint.set(ind);
    timePoint.time = R.timePoint.time(ind);

    timeInt.set = R.timeInterval.set(ind);
    timeInt.time = R.timeInterval.time(ind);
    if isfield(timeInt,'algebraic')
        timeInt.algebraic = R.timeInterval.algebraic(ind);
    end

    R = reachSet(timePoint,timeInt);
end

%------------- END OF CODE --------------