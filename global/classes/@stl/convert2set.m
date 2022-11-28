function res = convert2set(obj)
% convert2set - convert a logic formula to a set
%
% Syntax:  
%    res = convert2set(obj)
%
% Inputs:
%    obj - logic formula (class stl)
%
% Outputs:
%    res - resulting set (class halfspace or mptPolytope)
%
% Example: 
%    x = stl('x',2)
%    eq = x(1) < 5 & x(2) < 3;
%    set = convert2set(eq);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: stl

% Author:       Niklas Kochdumper
% Written:      9-November-2022 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    vars = getVariables(obj);
    res = convertToSetRecursive(obj,length(vars),vars);
end


% Auxiliary Functions -----------------------------------------------------

function res = convertToSetRecursive(obj,n,vars)
% recursive function to convert single parts of an STL-formula to a set

    switch obj.type

        case '&'
            set1 = convertToSetRecursive(obj.lhs,n,vars);
            set2 = convertToSetRecursive(obj.rhs,n,vars);
            res = mptPolytope(set1) & mptPolytope(set2);

        case '<'
            [c,d] = getHalfspace(obj.lhs,n,vars);
            res = halfspace(c,obj.rhs-d);

        case '<='
            [c,d] = getHalfspace(obj.lhs,n,vars);
            res = halfspace(c,obj.rhs-d);

        case '>'
            [c,d] = getHalfspace(obj.lhs,n,vars);
            res = halfspace(-c,-(obj.rhs-d));

        case '>='
            [c,d] = getHalfspace(obj.lhs,n,vars);
            res = halfspace(-c,-(obj.rhs-d));

        otherwise
            throw(CORAerror('CORA:notSupported',...
                  'Conversion to set not supported for these type of stl objects!'));
    end
end

function [c,d] = getHalfspace(obj,n,vars)
% convert a logic formula to a halfspace

    if isnumeric(obj)
        c = zeros(n,1); d = obj;
        return;
    end

    switch obj.type

        case 'variable'
            c = zeros(n,1); d = 0;
            ind = find(contains(vars,obj.var));
            c(ind(1)) = 1;

        case '+'
            [c1,d1] = getHalfspace(obj.lhs,n,vars);
            [c2,d2] = getHalfspace(obj.rhs,n,vars);
            c = c1 + c2; d = d1 + d2;

        case '-'
            [c1,d1] = getHalfspace(obj.lhs,n,vars);
            [c2,d2] = getHalfspace(obj.rhs,n,vars);
            c = c1 - c2; d = d1 - d2;

        case '*'
            [c,d] = getHalfspace(obj.rhs,n,vars);
            c = obj.lhs * c; d = obj.lhs*d;

        otherwise
            throw(CORAerror('CORA:notSupported',...
                  'Conversion to set not supported for these type of stl objects!'));
    end
end

%------------- END OF CODE --------------