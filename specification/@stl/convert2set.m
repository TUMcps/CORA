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
%    res - resulting set (class halfspace, polytope, or levelSet)
%
% Example: 
%    x = stl('x',2);
%    eq = x(1) < 5 & x(2) < 3;
%    set = convert2set(eq);
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

    vars = getVariables(obj);

    if isLinearPredicate(obj)
        res = aux_convertToSetRecursiveLinear(obj,length(vars),vars);
    else
        var = sym('x',[length(vars),1]);
        res = aux_convertToSetRecursiveNonlinear(obj,var,vars);
    end
end


% Auxiliary functions -----------------------------------------------------

function res = aux_convertToSetRecursiveLinear(obj,n,vars)
% recursive function to convert single parts of an STL-formula to a set

    switch obj.type

        case '&'
            set1 = aux_convertToSetRecursiveLinear(obj.lhs,n,vars);
            set2 = aux_convertToSetRecursiveLinear(obj.rhs,n,vars);
            res = and_(polytope(set1),polytope(set2),'exact');

        case '<'
            [c,d] = aux_getHalfspace(obj.lhs,n,vars);
            res = halfspace(c,obj.rhs-d);

        case '<='
            [c,d] = aux_getHalfspace(obj.lhs,n,vars);
            res = halfspace(c,obj.rhs-d);

        case '>'
            [c,d] = aux_getHalfspace(obj.lhs,n,vars);
            res = halfspace(-c,-(obj.rhs-d));

        case '>='
            [c,d] = aux_getHalfspace(obj.lhs,n,vars);
            res = halfspace(-c,-(obj.rhs-d));

        otherwise
            throw(CORAerror('CORA:notSupported',...
                  'Conversion to set not supported for these type of stl objects!'));
    end
end

function res = aux_convertToSetRecursiveNonlinear(obj,var,vars)
% recursive function to convert single parts of an STL-formula to a set

    switch obj.type

        case '&'
            set1 = aux_convertToSetRecursiveNonlinear(obj.lhs,var,vars);
            set2 = aux_convertToSetRecursiveNonlinear(obj.rhs,var,vars);
            res = set1 & set2;

        case '<'
            eq = aux_getLevelSet(obj.lhs,var,vars);
            res = levelSet(eq - obj.rhs,var,'<');

        case '<='
            eq = aux_getLevelSet(obj.lhs,var,vars);
            res = levelSet(eq - obj.rhs,var,'<=');

        case '>'
            eq = aux_getLevelSet(obj.lhs,var,vars);
            res = levelSet(-(eq - obj.rhs),var,'<');

        case '>='
            eq = aux_getLevelSet(obj.lhs,var,vars);
            res = levelSet(-(eq - obj.rhs),var,'<=');

        otherwise
            throw(CORAerror('CORA:notSupported',...
                  'Conversion to set not supported for these type of stl objects!'));
    end
end

function [c,d] = aux_getHalfspace(obj,n,vars)
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
            [c1,d1] = aux_getHalfspace(obj.lhs,n,vars);
            [c2,d2] = aux_getHalfspace(obj.rhs,n,vars);
            c = c1 + c2; d = d1 + d2;

        case '-'
            [c1,d1] = aux_getHalfspace(obj.lhs,n,vars);
            [c2,d2] = aux_getHalfspace(obj.rhs,n,vars);
            c = c1 - c2; d = d1 - d2;

        case '*'
            [c,d] = aux_getHalfspace(obj.rhs,n,vars);
            c = obj.lhs * c; d = obj.lhs*d;

        otherwise
            throw(CORAerror('CORA:notSupported',...
                  'Conversion to set not supported for these type of stl objects!'));
    end
end

function eq = aux_getLevelSet(obj,var,vars)
% convert a logic formula to a level set

    if isnumeric(obj)
        eq = obj;
        return;
    end

    switch obj.type

        case 'variable'
            ind = find(contains(vars,obj.var));
            eq = var(ind(1));

        case '+'
            eq1 = aux_getLevelSet(obj.lhs,var,vars);
            eq2 = aux_getLevelSet(obj.rhs,var,vars);
            eq = eq1 + eq2;

        case '-'
            eq1 = aux_getLevelSet(obj.lhs,var,vars);
            eq2 = aux_getLevelSet(obj.rhs,var,vars);
            eq = eq1 - eq2;

        case '*'
            eq1 = aux_getLevelSet(obj.lhs,var,vars);
            eq2 = aux_getLevelSet(obj.rhs,var,vars);
            eq = eq1 * eq2;
        
        case '^'
            eq1 = aux_getLevelSet(obj.lhs,var,vars);
            eq = eq1^obj.rhs;

        otherwise
            eq1 = aux_getLevelSet(obj.lhs,var,vars);
            eval(['eq = ',obj.type,'(eq1);']);
    end
end

% ------------------------------ END OF CODE ------------------------------
