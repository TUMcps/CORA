function [res,list] = assignIdentifiers(obj)
% assignIdentifiers - assign unique identifiers to non-temporal formulas
%
% Syntax:
%    res = assignIdentifiers(obj)
%
% Inputs:
%    obj - logic formula (class stl)
%
% Outputs:
%    res - resulting logic formula (class stl)
%    list - list of all non-temporal subformulas ordered by identifiers
%
% Example: 
%    x = stl('x',2);
%    eq = x(1) < 5 | globally(x(2) < 3,interval(0.1,0.2));
%    [eq_,list] = assignIdentifiers(eq)
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
    
    [res,~,list] = aux_recursive(obj,1,{});
end


% Auxiliary functions -----------------------------------------------------

function [res,cnt,list] = aux_recursive(obj,cnt,list)
% recursive function to assign ids to all non-temporal formulas

    if ~obj.temporal

        res = obj;
        res.id = cnt;
        list{cnt} = res;
        cnt = cnt + 1;

    elseif strcmp(obj.type,'next')
        
        [inner,cnt,list] = aux_recursive(obj.lhs,cnt,list);
        res = next(inner,obj.from);

    elseif strcmp(obj.type,'globally')

        [inner,cnt,list] = aux_recursive(obj.lhs,cnt,list);
        res = globally(inner,interval(obj.from,obj.to));

    elseif strcmp(obj.type,'finally')

        [inner,cnt,list] = aux_recursive(obj.lhs,cnt,list);
        res = finally(inner,interval(obj.from,obj.to));

    elseif strcmp(obj.type,'until')

        [lhs,cnt,list] = aux_recursive(obj.lhs,cnt,list);
        [rhs,cnt,list] = aux_recursive(obj.rhs,cnt,list);
        res = until(lhs,rhs,interval(obj.from,obj.to));

    elseif strcmp(obj.type,'release')

        [lhs,cnt,list] = aux_recursive(obj.lhs,cnt,list);
        [rhs,cnt,list] = aux_recursive(obj.rhs,cnt,list);
        res = release(lhs,rhs,interval(obj.from,obj.to));

    elseif strcmp(obj.type,'&')

        [lhs,cnt,list] = aux_recursive(obj.lhs,cnt,list);
        [rhs,cnt,list] = aux_recursive(obj.rhs,cnt,list);
        res = lhs & rhs;

    elseif strcmp(obj.type,'|')

        [lhs,cnt,list] = aux_recursive(obj.lhs,cnt,list);
        [rhs,cnt,list] = aux_recursive(obj.rhs,cnt,list);
        res = lhs | rhs;

    elseif strcmp(obj.type,'~')

        [inner,cnt,list] = aux_recursive(obj.lhs,cnt,list);
        res = ~inner;
    end  
end

% ------------------------------ END OF CODE ------------------------------
