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
% Last update:   08-February-2024 (FL, use interval property of stl instead of from and to)
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
        return
    end

    % temporal
    switch obj.type
        case 'next'
            % only compute lhs
            [inner,cnt,list] = aux_recursive(obj.lhs,cnt,list);
            res = next(inner,obj.from);

        case 'globally'
            % only compute lhs
            [inner,cnt,list] = aux_recursive(obj.lhs,cnt,list);
            res = globally(inner,obj.interval);
    
        case 'finally'
            % only compute lhs
            [inner,cnt,list] = aux_recursive(obj.lhs,cnt,list);
            res = finally(inner,obj.interval);
    
        case 'until'
            % compute both hs
            [lhs,cnt,list] = aux_recursive(obj.lhs,cnt,list);
            [rhs,cnt,list] = aux_recursive(obj.rhs,cnt,list);
            res = until(lhs,rhs,obj.interval);
    
        case 'release'
            % compute both hs
            [lhs,cnt,list] = aux_recursive(obj.lhs,cnt,list);
            [rhs,cnt,list] = aux_recursive(obj.rhs,cnt,list);
            res = release(lhs,rhs,obj.interval);
    
            % boolean
        case '&'
            % compute both hs
            [lhs,cnt,list] = aux_recursive(obj.lhs,cnt,list);
            [rhs,cnt,list] = aux_recursive(obj.rhs,cnt,list);
            res = lhs & rhs;
    
        case '|'
            % compute both hs
            [lhs,cnt,list] = aux_recursive(obj.lhs,cnt,list);
            [rhs,cnt,list] = aux_recursive(obj.rhs,cnt,list);
            res = lhs | rhs;
    
        case '~'
            % only compute lhs
            [inner,cnt,list] = aux_recursive(obj.lhs,cnt,list);
            res = ~inner;
    end  
end

% ------------------------------ END OF CODE ------------------------------
