function S = compact_(ls,method,tol,varargin)
% compact_ - aims at reducing the representation of a level set
%
% Syntax:
%    ls = compact_(ls)
%    ls = compact_(ls,method)
%    ls = compact_(ls,method,tol)
%
% Inputs:
%    ls - levelSet object
%    method - method for redundancy removal
%             'all' (default): checks for identical functions
%    tol - tolerance
%
% Outputs:
%    S - simplified representation (levelSet or emptySet)
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
% See also: contSet/compact

% Authors:       Mark Wetzlinger
% Written:       25-April-2023
% Last update:   ---
% Last revision: 30-July-2023 (MW, rename 'compact_')

% ------------------------------ BEGIN CODE -------------------------------

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

% ------------------------------ END OF CODE ------------------------------
