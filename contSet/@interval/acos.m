function I = acos(I)
% acos - Overloaded 'acos()' operator for intervals
%
% x_ is x infimum, x-- is x supremum
%
% [NaN, NaN] if (x_ < -1) and (x-- > 1),
% [NaN, NaN] if (x_ > 1) or (x-- < -1),
% [NaN, pi] if (x_ < -1) and (x-- in [-1, 1]),
% [0, NaN] if (x_ in [-1, 1]) and (x-- > 1),
% [acos(x--), acos(x_)] if (x >= -1) and (x <= 1).
%
% Syntax:  
%    I = acos(I)
%
% Inputs:
%    I - interval object
%
% Outputs:
%    res - interval object
%
% Example: 
%    I = interval([-0.5;0.3]);
%    res = acos(I);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Author:       Matthias Althoff
% Written:      05-February-2016
% Last update:  06-February-2016 (DG, Matrix case and typos)
%               20-February-2016 (DG, Errors are fixed, the matrix case is rewritten)
%               05-May-2020 (MW, standardized error message)
%               21-May-2022 (MW, remove new instantiation)
% Last revision:---

%------------- BEGIN CODE --------------

% scalar case
if isnumeric(I)

    if ((I.inf < -1) && (I.sup > 1)) || (I.inf > 1) || (I.sup < -1)
        I.inf = NaN;
        I.sup = NaN;
    elseif (I.inf < -1) && (I.sup < 1)
        I.inf = NaN;
        I.sup = pi;
    elseif (I.inf > -1) && (I.sup > 1)
        I.inf = 0;
        I.sup = NaN;
    else
        I.inf = acos(I.sup);
        I.sup = acos(I.inf);
    end
        
else

    lb = I.inf;
    ub = I.sup;
    
    % find indices
    ind1 = find(lb < -1 & ub > 1 | (lb > 1) | (ub < -1));
    I.inf(ind1) = NaN;
    I.sup(ind1) = NaN;
    
    ind2 = find(lb < -1 & ub >= -1 & ub <= 1);
    I.inf(ind2) = NaN;
    I.sup(ind2) = pi;
    
    ind3 = find(lb >= -1 & lb <= 1 & ub > 1);
    I.inf(ind3) = 0;
    I.sup(ind3) = NaN;
    
    ind4 = find(lb >= -1 & ub <= 1);
    I.inf(ind4) = acos(ub(ind4));
    I.sup(ind4) = acos(lb(ind4));
    
end

% return error if NaN occurs
if any(any(isnan(I.inf))) || any(any(isnan(I.sup)))
    throw(CORAerror('CORA:outOfDomain','validDomain','>= -1 && <= 1'));
end

%------------- END OF CODE --------------