function res = polytope2stl(obj,P)
% polytope2stl - convert polytope constraints to an STL-formula
%
% Syntax:
%    res = polytope2stl(obj,P)
%
% Inputs:
%    obj - logic formula (class stl)
%    P - polytope
%
% Outputs:
%    res - resulting stl formula (class stl)
%
% Example: 
%    x = stl('x',2)
%    P = polytope([1 1; -2 1; 0 -1],[1;1;1]);
%    res = polytope2stl(x,P);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: stl/in

% Authors:       Niklas Kochdumper, Benedikt Seidl
% Written:       09-November-2022 
% Last update:   ---
% Last revision: 23-September-2024 (MW, moved from global functions)

% ------------------------------ BEGIN CODE -------------------------------

    res = [];
    for i = 1:length(P.b)
        res = res & aux_halfspace2set(obj,P.A(i,:),P.b(i));
    end

end


% Auxiliary functions -----------------------------------------------------

function res = aux_halfspace2set(x,A,b)
% convert halfspace containment to a STL-formula
% TODO: simplify subsref syntax...

    ind = find(A ~= 0,1);
    if A(ind) == 1
        rhs = x.subsref(struct('type','()','subs',{{ind}}));
    else
        rhs = A(ind) * x.subsref(struct('type','()','subs',{{ind}}));
    end

    for i = ind+1:length(A)
        if A(i) ~= 0
            if A(i) == 1
                rhs = rhs + x.subsref(struct('type','()','subs',{{i}}));
            else
                rhs = rhs + A(i) * x.subsref(struct('type','()','subs',{{i}}));
            end
        end
    end

    res = rhs <= b;
end

% ------------------------------ END OF CODE ------------------------------
