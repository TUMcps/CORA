function I = log(I)
% log - Overloaded log (natural logorithm) function for intervals 
%
% x_ is x infimum, x-- is x supremum
%
% [NaN, NaN] if (x-- < 0),
% [NaN, log(x--)] if (x_ < 0) and (x-- >= 0)
% [log(x_), log(x--)] if (x_ >= 0).
%
% Syntax:
%    I = log(I)
%
% Inputs:
%    I - interval object
%
% Outputs:
%    I - interval object
%
% Example: 
%    I = interval([3;9]);
%    res = log(I)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Authors:       Dmitry Grebenyuk
% Written:       07-February-2016
% Last update:   21-February-2016 (DG, the matrix case is rewritten)
%                05-May-2020 (MW, addition of error message)
%                21-May-2022 (MW, remove new instantiation)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% scalar case
if isnumeric(I)
    
    if I.inf >= 0 
        I.inf = log(I.inf);
        I.sup = log(I.sup);
    elseif I.inf < 0 &&  I.sup >= 0
        I.inf = NaN;
        I.sup = log(I.sup);
    else
        I.inf = NaN;
        I.sup = NaN;
    end

else

    lb = I.inf;
    ub = I.sup;
    
    % find indices
    ind1 = find(lb >= 0);   
    I.inf(ind1) = log(lb(ind1));
    I.sup(ind1) = log(ub(ind1));
    
    ind2 = find(lb < 0 & ub >= 0);    
    I.inf(ind2) = NaN;
    I.sup(ind2) = log(ub(ind2));
    
    ind3 = find(ub < 0);    
    I.inf(ind3) = NaN;
    I.sup(ind3) = NaN;
       
end

% return error if NaN occures
if any(any(isnan(I.inf))) || any(any(isnan(I.sup)))
    throw(CORAerror('CORA:outOfDomain','validDomain','>= 0'));
end

% ------------------------------ END OF CODE ------------------------------
