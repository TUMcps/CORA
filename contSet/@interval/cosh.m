function I = cosh(I)
% cosh - Overloaded 'cosh()' operator for intervals
%
% x_ is x infimum, x-- is x supremum
%
% [cosh(x--), cosh(x_)] if (x-- <= 0),
% [cosh(x_), cosh(x--)] if (x_ >= 0),
% [1, max( acosh(x_), acosh(x--)] if (x_ < 0) and (x-- > 0).
%
% Syntax:
%    I = cosh(I)
%
% Inputs:
%    I - interval object
%
% Outputs:
%    I - interval object
%
% Example:
%    I = interval([-2;3],[3;4]);
%    res = cosh(I);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Authors:       Matthias Althoff
% Written:       05-February-2016
% Last update:   21-February-2016 (DG, the matrix case is rewritten)
%                21-May-2022 (MW, remove new instantiation)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% scalar case
if isnumeric(I)

    if I.sup <= 0
        I.inf = cosh(I.sup);
        I.sup = cosh(I.inf);
    elseif I.inf >= 0
        I.inf = cosh(I.inf);
        I.sup = cosh(I.sup);
    else
        I.inf = 1;
        I.sup = cosh( max( abs(I.inf), abs(I.sup) ) );
    end
        
else

    lb = I.inf;
    ub = I.sup;
    
    % find indices
    
    ind1 = find(ub <= 0);   
    I.inf(ind1) = cosh(ub(ind1));
    I.sup(ind1) = cosh(lb(ind1));
    
    ind2 = find(lb >= 0);    
    I.inf(ind2) = cosh(lb(ind2));
    I.sup(ind2) = cosh(ub(ind2));
    
    ind3 = find(lb <0 & ub > 0);    
    I.inf(ind3) = 1;
    I.sup(ind3) = cosh( max( abs(lb(ind3)), abs(ub(ind3)) ) );
       
end


% ------------------------------ END OF CODE ------------------------------
