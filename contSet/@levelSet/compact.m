function S = compact(ls)
% compact - aims at reducing the representation of a level set
%
% Syntax:  
%    ls = compact(ls)
%
% Inputs:
%    ls - levelSet object
%
% Outputs:
%    S - simplified representation
%
% Example: 
%    syms x y
%    eq1 = x^2 + y^2 - 4;
%    eq2 = -x^2 - y^2 + 4;
%    ls = levelSet([eq1;eq2],[x;y],{'<=','<='});
%
%    ls_ = compact(ls);   
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Mark Wetzlinger
% Written:      25-April-2023
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% assume no changes
S = ls;

if length(ls.eq) == 2
    
    % check if two equations are the same (-> reduce to only one)
    if logical(simplify(ls.eq(1) - ls.eq(2)) == 0)
        % check comparison operators
        if strcmp(ls.compOp{1},ls.compOp{2})
            % same comparison operators -> same equation
            S = levelSet(ls.eq(1),ls.vars,ls.compOp{1});        
    
        elseif ( strcmp(ls.compOp{1},'<=') || strcmp(ls.compOp{1},'==') ) ...
                && ( strcmp(ls.compOp{2},'<=') || strcmp(ls.compOp{2},'==') )
            % f(x) <= 0 and f(x) == 0  =>  f(x) == 0
            % note: prohibited by constructor!
            S = levelSet(ls.eq(1),ls.vars,'==');
    
        elseif ( strcmp(ls.compOp{1},'<=') || strcmp(ls.compOp{1},'<') ) ...
                && ( strcmp(ls.compOp{2},'<=') || strcmp(ls.compOp{2},'<') )
            % f(x) <= 0 and f(x) < 0  =>  f(x) < 0
            S = levelSet(ls.eq(1),ls.vars,'<');
    
        elseif ( strcmp(ls.compOp{1},'<') || strcmp(ls.compOp{1},'==') ) ...
                && ( strcmp(ls.compOp{2},'<') || strcmp(ls.compOp{2},'==') )
            % f(x) < 0 and f(x) == 0  =>  empty set
            % note: prohibited by constructor!
            S = emptySet(length(ls.vars));
        end
        return
    end
    
    % check if two equations are *(-1) of each other
    if logical(simplify(ls.eq(1) + ls.eq(2)) == 0)
    
        % check comparison operators ('<|<=' means either '<' or '<=' is used)
        if ( strcmp(ls.compOp{1},'<=') || strcmp(ls.compOp{1},'==') ) ...
                && ( strcmp(ls.compOp{2},'<=') || strcmp(ls.compOp{2},'==') )
            % f(x) <=|== 0 and -f(x) <=|== 0  =>  f(x) == 0
            S = levelSet(ls.eq(1),ls.vars,'==');
        
        elseif strcmp(ls.compOp{1},'<') || strcmp(ls.compOp{2},'<')
            % f(x) < 0 and -f(x) <=|<|== 0  
            %   or
            % f(x) <=|<|== 0 and -f(x) < 0  
            %   =>  empty set
            % note: some cases prohibited by constructor!
            S = emptySet(length(ls.vars));
        
        end
        % note: the above two if-branches should encompass all combinations
        return
    
    end

end

% other cases not yet investigated...

%------------- END OF CODE --------------