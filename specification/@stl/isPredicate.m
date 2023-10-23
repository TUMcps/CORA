function res = isPredicate(obj)
% isPredicate - check if the STL object represents a predicate
%
% Syntax:
%    res = isPredicate(obj)
%
% Inputs:
%    obj - logic formula (class stl)
%
% Outputs:
%    res - is predicate (res = true) or not (res = false)
%
% Example:
%    x = stl('x',2);
%    eq1 = x(1) < 5 & x(2) > 2;
%    eq2 = x(1)*x(2) < 3;
%    isPredicate(eq1)
%    isPredicate(eq2)
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

    elseif ismember(obj.type,{'&','<=','<','>','>=','+','-','*','^'})

        res = aux_recursive(obj.lhs) & aux_recursive(obj.rhs);

    else

        res = false;

    end
end

% ------------------------------ END OF CODE ------------------------------
