function res = minimumTime(obj)
% minimumTime - get minimum time for evaluation the temporal logic formula
%
% Syntax:
%    res = minimumTime(obj)
%
% Inputs:
%    obj - temporal logic formula (class stl)
%
% Outputs:
%    res - minimum time required for evaluating the temporal logic formula
% 
% Example: 
%    x = stl('x',2);
%    eq = next(globally(finally(x(2) < 3,interval(0,7)),interval(1,3)),2);
%    time = minimumTime(eq)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: stl

% Authors:       Niklas Kochdumper, Benedikt Seidl
% Written:       09-November-2022 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------
    
   res = aux_recursive(obj,0);

end


% Auxiliary functions -----------------------------------------------------

function res = aux_recursive(obj,time)
% recursive function for computing the maximum time

    if ~obj.temporal

        res = time;

    elseif strcmp(obj.type,'next')

        res = aux_recursive(obj.lhs,obj.from + time);

    elseif strcmp(obj.type,'~') 

        res = aux_recursive(obj.lhs,time);

    elseif ismember(obj.type,{'&','|'})

        res = min(aux_recursive(obj.lhs,time),aux_recursive(obj.rhs,time));

    elseif ismember(obj.type,{'finally','globally'})

        res = obj.from + time;

    elseif ismember(obj.type,{'until','release'})

        res = time;

    end
end

% ------------------------------ END OF CODE ------------------------------
