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

% Authors:       Niklas Kochdumper
% Written:       02-June-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check input arguments (only first two)
inputArgsCheck({{R,'att','reachSet'};
                {prop,'str',{'location','parent','time'}}});

% init result
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

    res = [];

    if ~isa(val,'interval')
        val = interval(val-1e-10,val+1e-10);
    end

    for i = 1:size(R,1)

       timeInt = []; timePoint = [];

       if ~isempty(R(i).timeInterval)
       
           ind = [];
           
           for j = 1:length(R(i).timeInterval.time) 
              if isIntersecting_(R(i).timeInterval.time{j},val,'exact')
                  ind = [ind;j];
              end
           end
    
           if ~isempty(ind)
               timeInt.set = R(i).timeInterval.set(ind);
               timeInt.time = R(i).timeInterval.time(ind);
               if isfield(R(i).timeInterval,'algebraic')
                   timeInt.algebraic = R.timeInterval.algebraic(ind);
               end
           end
       end

       ind = [];
           
       for j = 1:length(R(i).timePoint.time) 
          if isIntersecting_(interval(R(i).timePoint.time{j}),val,'exact')
              ind = [ind;j];
          end
       end

       if ~isempty(ind)
           timePoint.set = R(i).timePoint.set(ind);
           timePoint.time = R(i).timePoint.time(ind);
       end

       if ~isempty(timeInt)
           res = add(res,reachSet(timePoint,timeInt)); 
       elseif ~isempty(timePoint)
           res = add(res,reachSet(timePoint));
       end   
    end
end

end


% Auxiliary functions -----------------------------------------------------

function R = aux_getSubset(R,ind)
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

% ------------------------------ END OF CODE ------------------------------
