function res = acosh(intVal)
% acosh - Overloaded 'acosh()' operator for intervals
%
% Syntax:  
%    res = acosh(intVal)
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
% Last update:  21-February-2016 the matrix case is rewritten  (Dmitry Grebenyuk)
% Last revision:---

% x_ is x infimum, x-- is x supremum

% [NaN, NaN] if (x-- < 1),
% [NaN, acosh(x--)] if (x_ < 1) and (x-- >= 1),
% [acosh(x_), acosh(x--)] if (x_ >= 1).

%------------- BEGIN CODE --------------

% scalar case
if isnumeric(intVal)
    
        res = interval();

        if intVal.inf >= 1
            res.inf = acosh(intVal.inf);
            res.sup = acosh(intVal.sup);
        elseif intVal.sup >= 1
            res.inf = NaN;
            res.sup = acosh(intVal.sup);
        else
            res.inf = NaN;
            res.sup = NaN;
        end
   
else

    % to preserve the shape    
    res = intVal;
    
    % find indices
    
    ind1 = find(intVal.inf >= 1);   
    res.inf(ind1) = acosh(intVal.inf(ind1));
    res.sup(ind1) = acosh(intVal.sup(ind1));
    
    ind2 = find(intVal.inf < 1 & intVal.sup >= 1);    
    res.inf(ind2) = NaN;
    res.sup(ind2) = acosh(intVal.sup(ind2));
    
    ind3 = find(intVal.sup < 1);    
    res.inf(ind3) = NaN;
    res.sup(ind3) = NaN;
       
end

%------------- END OF CODE --------------