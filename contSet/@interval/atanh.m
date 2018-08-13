function res = atanh(intVal)
% atanh - Overloaded 'atanh()' operator for intervals
%
% Syntax:  
%    res = atanh(intVal)
%
% Inputs:
%    intVal - interval object
%
% Outputs:
%    res - interval object
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Author:       Matthias Althoff
% Written:      12-February-2016
% Last update:  21-February-2016 Errors are fixed, the matrix case is rewritten  (Dmitry Grebenyuk)
% Last revision:---

% x_ is x infimum, x-- is x supremum

% [NaN, NaN] if (x_ < -1) and (x-- > 1),
% [NaN, NaN] if (x_ > 1) or (x-- < -1),
% [NaN, atanh(x--)] if (x_ < -1) and (x-- in [-1, 1]),
% [atanh(x_), NaN] if (x_ in [-1, 1]) and (x-- > 1),
% [atanh(x_), atanh(x--)] if (x >= -1) and (x <= 1).

%------------- BEGIN CODE --------------

% scalar case
if isnumeric(intVal)
    
        res = interval();

        if ((intVal.inf < -1) && (intVal.sup > 1)) || (intVal.inf > 1) || (intVal.sup < -1)
            res.inf = NaN;
            res.sup = NaN;
        elseif (intVal.inf < -1) && (intVal.sup <= 1)
            res.inf = NaN;
            res.sup = atanh(intVal.sup);
        elseif (intVal.inf >= -1) && (intVal.sup > 1)
            res.inf = atanh(intVal.inf);
            res.sup = NaN;
        else
            res.inf = atanh(intVal.inf);
            res.sup = atanh(intVal.sup);
        end
        
else

% matrix case
% rand(100, 100)
%time_CORA =
%   0.001914723042715
%time_INTLAB =
%   0.052251330145747

    % to preserve the shape    
    res = intVal;
    
    % find indices
    
    ind1 = find((intVal.inf < -1 & intVal.sup > 1) | (intVal.inf > 1) | (intVal.sup < -1));   
    res.inf(ind1) = NaN;
    res.sup(ind1) = NaN;
    
    ind2 = find(intVal.inf < -1 & intVal.sup >= -1 & intVal.sup <= 1);    
    res.inf(ind2) = NaN;
    res.sup(ind2) = atanh(intVal.sup(ind2));
    
    ind3 = find(intVal.inf >= -1 & intVal.inf <= 1 & intVal.sup > 1);    
    res.inf(ind3) = atanh(intVal.inf(ind3));
    res.sup(ind3) = NaN;
    
    ind4 = find(intVal.inf >= -1 & intVal.sup <= 1);    
    
    res.inf(ind4) = atanh(intVal.inf(ind4));
    res.sup(ind4) = atanh(intVal.sup(ind4));
       

end
%------------- END OF CODE --------------