function res = isLinearPredicate(obj)
% isLinearPredicate - check if the STL object represents a linear predicate
%
% Syntax:
%    res = isLinearPredicate(obj)
%
% Inputs:
%    obj - logic formula (class stl)
%
% Outputs:
%    res - is linear predicate (res = true) or not (res = false)
%
% Example: 
%    x = stl('x',2);
%    eq1 = x(1) < 5 & x(2) > 2;
%    eq2 = x(1)*x(2) < 3;
%    isLinearPredicate(eq1)
%    isLinearPredicate(eq2)
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
    
    if obj.temporal
        res = false; return;
    end

    res = aux_recursive(obj);
end


% Auxiliary functions -----------------------------------------------------

function res = aux_recursive(obj)
% recursive function to check object is a linear predicate

    if isnumeric(obj) || strcmp(obj.type,'variable')

        res = true;

    elseif ismember(obj.type,{'&','<=','<','>','>=','+','-'})

        res = aux_recursive(obj.lhs) & aux_recursive(obj.rhs);

    elseif strcmp(obj.type,'*')

        if isnumeric(obj.lhs)
            res = aux_recursive(obj.rhs);
        elseif isnumeric(obj.rhs)
            res = aux_recursive(obj.lhs);
        else
            res = false;
        end

    else

        res = false;

    end
end

% ------------------------------ END OF CODE ------------------------------
