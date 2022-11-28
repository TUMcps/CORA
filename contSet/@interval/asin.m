function I = asin(I)
% asin - Overloaded 'asin()' operator for intervals
%
% x_ is x infimum, x-- is x supremum
%
% [NaN, NaN] if (x_ < -1) and (x-- > 1),
% [NaN, NaN] if (x_ > 1) or (x-- < -1),
% [NaN, asin(x--)] if (x_ < -1) and (x-- in [-1, 1]),
% [asin(x_), NaN] if (x_ in [-1, 1]) and (x-- > 1),
% [asin(x_), asin(x--)] if (x >= -1) and (x <= 1).
%
% Syntax:  
%    I = asin(I)
%
% Inputs:
%    I - interval object
%
% Outputs:
%    I - interval object
%
% Example:
%    I = interval([-0.5;0.3]);
%    res = asin(I);
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
        I.sup = asin(I.sup);
    elseif (I.inf > -1) && (I.sup > 1)
        I.inf = asin(I.inf);
        I.sup = NaN;
    else
        I.inf = asin(I.inf);
        I.sup = asin(I.sup);
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
    I.sup(ind2) = asin(ub(ind2));
    
    ind3 = find(lb >= -1 & lb <= 1 & ub > 1);    
    I.inf(ind3) = asin(lb(ind3));
    I.sup(ind3) = NaN;
    
    ind4 = find(lb >= -1 & ub <= 1);    
    
    I.inf(ind4) = asin(lb(ind4));
    I.sup(ind4) = asin(ub(ind4));

end

% return error if NaN occures
if any(any(isnan(I.inf))) || any(any(isnan(I.sup)))
    throw(CORAerror('CORA:outOfDomain','validDomain','>= -1 && <= 1'));
end


%------------- END OF CODE --------------