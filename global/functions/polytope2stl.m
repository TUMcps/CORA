function res = polytope2stl(obj,set)
% polytope - convert polytope containment constraint to a STL-formula
%
% Syntax:  
%    res = polytope2stl(obj,set)
%
% Inputs:
%    obj - logic formula (class stl)
%    set - set (class mptPolytope)
%
% Outputs:
%    res - resulting stl formula (class stl)
%
% Example: 
%    x = stl('x',2)
%    pgon = polygon.generateRandom();
%    res = polytope2stl(x,set
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: stl/in

% Author:       Niklas Kochdumper
% Written:      9-November-2022 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    res = [];

    for i = 1:size(set.P.A,1)
        res = res & halfspace2set(obj,set.P.A(i,:),set.P.b(i));
    end
end


% Auxiliary Functions -----------------------------------------------------

function res = halfspace2set(x,c,d)
% convert halfspace containment to a STL-formula

    rhs = c(1) * x(1);

    for i = 2:length(c)
        rhs = rhs + c(i) * x(i);
    end

    res = rhs <= d;
end

%------------- END OF CODE --------------