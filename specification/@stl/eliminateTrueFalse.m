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
%    x = stl('x',2);
%    eq = until(x(1) < 5 | true,x(2) > 5,interval(0.2,0.3));
%    eq_ = eliminateTrueFalse(eq)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: stl

% Authors:       Niklas Kochdumper, Benedikt Seidl
% Written:       09-November-2022 
% Last update:   08-February-2024 (FL, use interval property of stl instead of from and to)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    if ~isempty(obj.id)

        res = obj;

        % boolean
    elseif strcmp(obj.type,'&') % ---

        lhs = eliminateTrueFalse(obj.lhs);
        rhs = eliminateTrueFalse(obj.rhs);

        if strcmp(lhs.type,'false') || strcmp(rhs.type,'false')
            res = stl(false);
        elseif strcmp(lhs.type,'true') 
            res = rhs;
        elseif strcmp(rhs.type,'true')
            res = lhs;
        else
            res = lhs & rhs;
        end

    elseif strcmp(obj.type,'|') % ---

        lhs = eliminateTrueFalse(obj.lhs);
        rhs = eliminateTrueFalse(obj.rhs);

        if strcmp(lhs.type,'true') || strcmp(rhs.type,'true')
            res = stl(true);
        elseif strcmp(lhs.type,'false') 
            res = rhs;
        elseif strcmp(rhs.type,'false')
            res = lhs;
        else
            res = lhs | rhs;
        end

    elseif strcmp(obj.type,'~') % ---

        inner = eliminateTrueFalse(obj.lhs);

        if strcmp(inner.type,'false')
            res = stl(true);
        elseif strcmp(inner.type,'true')
            res = stl(false);
        else
            res = ~inner;
        end

        % temporal
    elseif strcmp(obj.type,'next') % ---

        inner = eliminateTrueFalse(obj.lhs);

        if strcmp(inner.type,'false')
            res = stl(false);
        elseif strcmp(inner.type,'true')
            res = stl(true);
        else
            res = next(inner,obj.from);
        end

    elseif strcmp(obj.type,'until') % ---

        lhs = eliminateTrueFalse(obj.lhs);
        rhs = eliminateTrueFalse(obj.rhs);

        if strcmp(lhs.type,'true')
            res = finally(rhs,obj.interval);
        elseif strcmp(rhs.type,'true')
            res = stl(true);
        elseif strcmp(rhs.type,'false')
            res = stl(false);
        else
            res = until(lhs,rhs,obj.interval);
        end

    elseif strcmp(obj.type,'release') % ---

        lhs = eliminateTrueFalse(obj.lhs);
        rhs = eliminateTrueFalse(obj.rhs);

        if strcmp(lhs.type,'false')
            res = globally(rhs,obj.interval);
        elseif strcmp(rhs.type,'true')
            res = stl(true);
        elseif strcmp(lhs.type,'true')
            res = stl(true);
        else
            res = release(lhs,rhs,obj.interval);
        end

    elseif strcmp(obj.type,'finally') % ---

        inner = eliminateTrueFalse(obj.lhs);
        
        if strcmp(inner.type,'false')
            res = stl(false);
        elseif strcmp(inner.type,'true')
            res = stl(true);
        else
            res = finally(inner,obj.interval);
        end

    elseif strcmp(obj.type,'globally') % ---

        inner = eliminateTrueFalse(obj.lhs);
        
        if strcmp(inner.type,'false')
            res = stl(false);
        elseif strcmp(inner.type,'true')
            res = stl(true);
        else
            res = globally(inner,obj.interval);
        end    

        % comparator
    elseif ismember(obj.type,{'<','<=','>','>=','true','false'})
        res = obj;
    end
end

% ------------------------------ END OF CODE ------------------------------
