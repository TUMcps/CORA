function res = cosh(intVal)
% cosh - Overloaded 'cosh()' operator for intervals
%
% Syntax:  
%    res = cosh(intVal)
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
% Written:      05-February-2016
% Last update:  21-February-2016 the matrix case is rewritten  (Dmitry Grebenyuk)
% Last revision:---

% x_ is x infimum, x-- is x supremum

% [cosh(x--), cosh(x_)] if (x-- <= 0),
% [cosh(x_), cosh(x--)] if (x_ >= 0),
% [1, max( acosh(x_), acosh(x--)] if (x_ < 0) and (x-- > 0).

%------------- BEGIN CODE --------------

% scalar case
if isnumeric(intVal)

        res = interval();

        if intVal.sup <= 0
            res.inf = cosh(intVal.sup);
            res.sup = cosh(intVal.inf);
        elseif intVal.inf >= 0
            res.inf = cosh(intVal.inf);
            res.sup = cosh(intVal.sup);
        else
            res.inf = 1;
            res.sup = cosh( max( abs(intVal.inf), abs(intVal.sup) ) );
        end
        
else

% matrix case
% rand(100, 100)
%time_CORA =
%   0.001808342961487
%time_INTLAB =
%   0.006772680963571

    % to preserve the shape    
    res = intVal;
    
    % find indices
    
    ind1 = find(intVal.sup <= 0);   
    res.inf(ind1) = cosh(intVal.sup(ind1));
    res.sup(ind1) = cosh(intVal.inf(ind1));
    
    ind2 = find(intVal.inf >= 0);    
    res.inf(ind2) = cosh(intVal.inf(ind2));
    res.sup(ind2) = cosh(intVal.sup(ind2));
    
    ind3 = find(intVal.inf <0 & intVal.sup > 0);    
    res.inf(ind3) = 1;
    res.sup(ind3) = cosh( max( abs(intVal.inf(ind3)), abs(intVal.sup(ind3)) ) );
       
end


%------------- END OF CODE --------------