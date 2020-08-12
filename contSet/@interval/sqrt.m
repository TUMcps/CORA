function res = sqrt(intVal)
% sqrt - Overloaded sqrt function for intervals 
%
% x_ is x infimum, x-- is x supremum
%
% [NaN, NaN] if (x-- < 0),
% [NaN, sqrt(x--)] if (x_ < 0) and (x-- >= 0),
% [sqrt(x_), sqrt(x--)] if (x_ >= 0).
%
% Syntax:  
%    res = sqrt(intVal)
%
% Inputs:
%    intVal - argument of which square root should be obtained
%
% Outputs:
%    res - interval
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Author:       Matthias Althoff
% Written:      20-January-2016
%               21-February-2016 (DG, the matrix case is rewritten)
%               05-May-2020 (MW, standardized error message)
% Last revision:---

%------------- BEGIN CODE --------------

% scalar case
if isnumeric(intVal)

    res = interval();

    if intVal.inf >= 0 
        res.inf = sqrt(intVal.inf);
        res.sup = sqrt(intVal.sup);
    elseif intVal.inf < 0 &&  intVal.sup >= 0
        res.inf = NaN;
        res.sup = sqrt(intVal.sup);
    else
        res.inf = NaN;
        res.sup = NaN;
    end

else

    % to preserve the shape    
    res = intVal;
    
    % find indices
    
    ind1 = find(intVal.inf >= 0);   
    res.inf(ind1) = sqrt(intVal.inf(ind1));
    res.sup(ind1) = sqrt(intVal.sup(ind1));
    
    ind2 = find(intVal.inf < 0 & intVal.sup >= 0);    
    res.inf(ind2) = NaN;
    res.sup(ind2) = sqrt(intVal.sup(ind2));
    
    ind3 = find(intVal.sup < 0);    
    res.inf(ind3) = NaN;
    res.sup(ind3) = NaN;      

end

% return error if NaN occures
if any(any(isnan(res.inf))) || any(any(isnan(res.sup)))
	error(resIntNaNInf()); 
end

%------------- END OF CODE --------------