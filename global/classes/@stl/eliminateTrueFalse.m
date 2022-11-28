function res = eliminateTrueFalse(obj)
% eliminateTrueFalse - remove true and false from STL formulas 
%
% Syntax:  
%    res = eliminateTrueFalse(obj)
%
% Inputs:
%    obj - logic formula (class stl)
%
% Outputs:
%    res - resulting logic formula (class stl)
%
% Example: 
%    x = stl('x',2)
%    eq = until(x(1) < 5 | true,x(2) > 5,interval(0.2,0.3));
%    eq_ = eliminateTrueFalse(eq)
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

    if strcmp(obj.type,'&')

        lhs = eliminateTrueFalse(obj.lhs);
        rhs = eliminateTrueFalse(obj.rhs);

        if strcmp(lhs.type,'false') || strcmp(rhs.type,'false')
            res = createBoolean(false);
        elseif strcmp(lhs.type,'true') 
            res = rhs;
        elseif strcmp(rhs.type,'true')
            res = lhs;
        else
            res = lhs & rhs;
        end

    elseif strcmp(obj.type,'|')

        lhs = eliminateTrueFalse(obj.lhs);
        rhs = eliminateTrueFalse(obj.rhs);

        if strcmp(lhs.type,'true') || strcmp(rhs.type,'true')
            res = createBoolean(true);
        elseif strcmp(lhs.type,'false') 
            res = rhs;
        elseif strcmp(rhs.type,'false')
            res = lhs;
        else
            res = lhs | rhs;
        end

    elseif strcmp(obj.type,'~')

        inner = eliminateTrueFalse(obj.lhs);

        if strcmp(inner.type,'false')
            res = createBoolean(true);
        elseif strcmp(inner.type,'true')
            res = createBoolean(false);
        else
            res = ~inner;
        end

    elseif strcmp(obj.type,'next')

        inner = eliminateTrueFalse(obj.lhs);

        if strcmp(inner.type,'false')
            res = createBoolean(false);
        elseif strcmp(inner.type,'true')
            res = createBoolean(true);
        else
            res = next(inner,obj.from);
        end

    elseif strcmp(obj.type,'until')

        lhs = eliminateTrueFalse(obj.lhs);
        rhs = eliminateTrueFalse(obj.rhs);

        if strcmp(lhs.type,'true')
            res = finally(rhs,interval(obj.from,obj.to));
        elseif strcmp(rhs.type,'true')
            res = createBoolean(true);
        elseif strcmp(rhs.type,'false')
            res = createBoolean(false);
        else
            res = until(lhs,rhs,interval(obj.from,obj.to));
        end

    elseif strcmp(obj.type,'release')

        lhs = eliminateTrueFalse(obj.lhs);
        rhs = eliminateTrueFalse(obj.rhs);

        if strcmp(lhs.type,'false')
            res = globally(rhs,interval(obj.from,obj.to));
        elseif strcmp(rhs.type,'true')
            res = createBoolean(true);
        elseif strcmp(lhs.type,'true')
            res = createBoolean(true);
        else
            res = release(lhs,rhs,interval(obj.from,obj.to));
        end

    elseif strcmp(obj.type,'finally')

        inner = eliminateTrueFalse(obj.lhs);
        
        if strcmp(inner.type,'false')
            res = createBoolean(false);
        elseif strcmp(inner.type,'true')
            res = createBoolean(true);
        else
            res = finally(inner,interval(obj.from,obj.to));
        end

    elseif strcmp(obj.type,'globally')

        inner = eliminateTrueFalse(obj.lhs);
        
        if strcmp(inner.type,'false')
            res = createBoolean(false);
        elseif strcmp(inner.type,'true')
            res = createBoolean(true);
        else
            res = globally(inner,interval(obj.from,obj.to));
        end    

    elseif ismember(obj.type,{'<','<=','>','>=','true','false'})
        res = obj;
    end
end


% Auxiliary Functions -----------------------------------------------------

function obj = createBoolean(val)
% convert a boolean value to a stl object

    if val
        obj = stl('x',1);
        obj.type = 'true';
    else
        obj = stl('x',1);
        obj.type = 'false';
    end
    obj.variables = [];
    obj.logic = true;
end

%------------- END OF CODE --------------