function res = negationNormalForm(obj)
% negationNormalForm - convert STL formula to negation normal form
%
% Syntax:  
%    res = negationNormalForm(obj)
%
% Inputs:
%    obj - logic formula (class stl)
%
% Outputs:
%    res - resulting stl formula in negation normal form (class stl)
%
% Example: 
%    x = stl('x',2)
%    eq = ~(x(1) < 5 | globally(x(2) < 3,interval(0.1,0.2)));
%    eq_ = negationNormalForm(eq)
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

    res = recursive(obj,false);
    
end

% Auxiliary Functions -----------------------------------------------------

function res = recursive(obj,neg)
% recursive function to convert single parts of an STL formula to negation
% normal form

    % negate formula (neg = true) or not (neg = false)
    if neg

        if strcmp(obj.type,'&')
            res = recursive(obj.lhs,true) | recursive(obj.rhs,true);
        elseif strcmp(obj.type,'|')
            res = recursive(obj.lhs,true) & recursive(obj.rhs,true);
        elseif strcmp(obj.type,'~')
            res = recursive(obj.lhs,false);
        elseif strcmp(obj.type,'finally')
            res = globally(recursive(obj.lhs,true), ...
                                            interval(obj.from,obj.to));
        elseif strcmp(obj.type,'globally')
            res = finally(recursive(obj.lhs,true), ...
                                            interval(obj.from,obj.to));
        elseif strcmp(obj.type,'until')
            res = release(recursive(obj.lhs,true), ...
                        recursive(obj.rhs,true),interval(obj.from,obj.to));
        elseif strcmp(obj.type,'release')
            res = until(recursive(obj.lhs,true), ...
                        recursive(obj.rhs,true),interval(obj.from,obj.to));
        elseif strcmp(obj.type,'next')
            res = next(recursive(obj.lhs,true),obj.from);
        elseif strcmp(obj.type,'<')
            res = obj.lhs > obj.rhs;
        elseif strcmp(obj.type,'>')
            res = obj.lhs < obj.rhs;
        elseif strcmp(obj.type,'<=')
            res = obj.lhs >= obj.rhs;
        elseif strcmp(obj.type,'>=')
            res = obj.lhs <= obj.rhs;
        end

    else

        if ~obj.temporal && ~ismember(obj.type,{'&','|','~'})
            res = obj;
        elseif strcmp(obj.type,'&')
            res = recursive(obj.lhs,false) & recursive(obj.rhs,false);
        elseif strcmp(obj.type,'|')
            res = recursive(obj.lhs,false) | recursive(obj.rhs,false);
        elseif strcmp(obj.type,'~')
            res = recursive(obj.lhs,true);
        elseif strcmp(obj.type,'finally')
            res = finally(recursive(obj.lhs,false), ...
                                            interval(obj.from,obj.to));
        elseif strcmp(obj.type,'globally')
            res = globally(recursive(obj.lhs,false), ...
                                            interval(obj.from,obj.to));
        elseif strcmp(obj.type,'until')
            res = until(recursive(obj.lhs,false), ...
                       recursive(obj.rhs,false),interval(obj.from,obj.to));
        elseif strcmp(obj.type,'release')
            res = release(recursive(obj.lhs,false), ...
                       recursive(obj.rhs,false),interval(obj.from,obj.to));
        elseif strcmp(obj.type,'next')
            res = next(recursive(obj.lhs,false),obj.from);
        end
    end
end

%------------- END OF CODE --------------