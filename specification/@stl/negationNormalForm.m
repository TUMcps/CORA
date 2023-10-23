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
%    x = stl('x',2);
%    eq = ~(x(1) < 5 | globally(x(2) < 3,interval(0.1,0.2)));
%    eq_ = negationNormalForm(eq)
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

    res = aux_recursive(obj,false);
    
end


% Auxiliary functions -----------------------------------------------------

function res = aux_recursive(obj,neg)
% recursive function to convert single parts of an STL formula to negation
% normal form

    % negate formula (neg = true) or not (neg = false)
    if neg

        if strcmp(obj.type,'&')
            res = aux_recursive(obj.lhs,true) | aux_recursive(obj.rhs,true);
        elseif strcmp(obj.type,'|')
            res = aux_recursive(obj.lhs,true) & aux_recursive(obj.rhs,true);
        elseif strcmp(obj.type,'~')
            res = aux_recursive(obj.lhs,false);
        elseif strcmp(obj.type,'finally')
            res = globally(aux_recursive(obj.lhs,true), ...
                                            interval(obj.from,obj.to));
        elseif strcmp(obj.type,'globally')
            res = finally(aux_recursive(obj.lhs,true), ...
                                            interval(obj.from,obj.to));
        elseif strcmp(obj.type,'until')
            res = release(aux_recursive(obj.lhs,true), ...
                        aux_recursive(obj.rhs,true),interval(obj.from,obj.to));
        elseif strcmp(obj.type,'release')
            res = until(aux_recursive(obj.lhs,true), ...
                        aux_recursive(obj.rhs,true),interval(obj.from,obj.to));
        elseif strcmp(obj.type,'next')
            res = next(aux_recursive(obj.lhs,true),obj.from);
        elseif strcmp(obj.type,'<')
            res = obj.lhs > obj.rhs;
        elseif strcmp(obj.type,'>')
            res = obj.lhs < obj.rhs;
        elseif strcmp(obj.type,'<=')
            res = obj.lhs >= obj.rhs;
        elseif strcmp(obj.type,'>=')
            res = obj.lhs <= obj.rhs;
        elseif strcmp(obj.type,'true')
            res = stl(false);
        elseif strcmp(obj.type,'false')
            res = stl(true);
        else
            res = ~obj;
        end

    else

        if strcmp(obj.type,'&')
            res = aux_recursive(obj.lhs,false) & aux_recursive(obj.rhs,false);
        elseif strcmp(obj.type,'|')
            res = aux_recursive(obj.lhs,false) | aux_recursive(obj.rhs,false);
        elseif strcmp(obj.type,'~')
            res = aux_recursive(obj.lhs,true);
        elseif strcmp(obj.type,'finally')
            res = finally(aux_recursive(obj.lhs,false), ...
                                            interval(obj.from,obj.to));
        elseif strcmp(obj.type,'globally')
            res = globally(aux_recursive(obj.lhs,false), ...
                                            interval(obj.from,obj.to));
        elseif strcmp(obj.type,'until')
            res = until(aux_recursive(obj.lhs,false), ...
                       aux_recursive(obj.rhs,false),interval(obj.from,obj.to));
        elseif strcmp(obj.type,'release')
            res = release(aux_recursive(obj.lhs,false), ...
                       aux_recursive(obj.rhs,false),interval(obj.from,obj.to));
        elseif strcmp(obj.type,'next')
            res = next(aux_recursive(obj.lhs,false),obj.from);
        else
            res = obj;
        end

    end
end

% ------------------------------ END OF CODE ------------------------------
