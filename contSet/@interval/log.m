function res = log(intVal)
% sqrt - Overloaded log (natural logorithm) function for intervals 
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

% Author:       Dmitry Grebenyuk
% Written:      07-February-2016
% Last update:  21-February-2016 the matrix case is rewritten  (Dmitry Grebenyuk)
% Last revision:---

% x_ is x infimum, x-- is x supremum

% [NaN, NaN] if (x-- < 0),
% [NaN, log(x--)] if (x_ < 0) and (x-- >= 0)
% [log(x_), log(x--)] if (x_ >= 0).
%------------- BEGIN CODE --------------

% scalar case
if isnumeric(intVal)
    
        res = interval();

        if intVal.inf >= 0 
            res.inf = log(intVal.inf);
            res.sup = log(intVal.sup);
        elseif intVal.inf < 0 &&  intVal.sup >= 0
            res.inf = NaN;
            res.sup = log(intVal.sup);
        else
            res.inf = NaN;
            res.sup = NaN;
        end
        
else

% matrix case
% rand(100, 100)
%time_CORA =
%   0.001682028950382
%time_INTLAB =
%   0.015084498152460

    % to preserve the shape    
    res = intVal;
    
    % find indices
    
    ind1 = find(intVal.inf >= 0);   
    res.inf(ind1) = log(intVal.inf(ind1));
    res.sup(ind1) = log(intVal.sup(ind1));
    
    ind2 = find(intVal.inf < 0 & intVal.sup >= 0);    
    res.inf(ind2) = NaN;
    res.sup(ind2) = log(intVal.sup(ind2));
    
    ind3 = find(intVal.sup < 0);    
    res.inf(ind3) = NaN;
    res.sup(ind3) = NaN;
       
end


%------------- END OF CODE --------------
