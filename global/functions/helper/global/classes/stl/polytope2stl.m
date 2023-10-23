function res = polytope2stl(obj,set)
% polytope2stl - convert polytope containment constraint to a STL-formula
%
% Syntax:
%    res = polytope2stl(obj,set)
%
% Inputs:
%    obj - logic formula (class stl)
%    set - set (class polytope)
%
% Outputs:
%    res - resulting stl formula (class stl)
%
% Example: 
%    x = stl('x',2)
%    set = polytope([1 1; -2 1; 0 -1],[1;1;1]);
%    res = polytope2stl(x,set);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: stl/in

% Authors:       Niklas Kochdumper, Benedikt Seidl
% Written:       09-November-2022 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    res = [];
    for i = 1:length(set.b)
        res = res & aux_halfspace2set(obj,set.A(i,:),set.b(i));
    end

end


% Auxiliary functions -----------------------------------------------------

function res = aux_halfspace2set(x,c,d)
% convert halfspace containment to a STL-formula

    ind = find(c ~= 0,1);
    if c(ind) == 1
        rhs = x(ind);
    else
        rhs = c(ind) * x(ind);
    end

    for i = ind+1:length(c)
        if c(i) ~= 0
            if c(i) == 1
                rhs = rhs + x(i);
            else
                rhs = rhs + c(i) * x(i);
            end
        end
    end

    res = rhs <= d;
end

% ------------------------------ END OF CODE ------------------------------
