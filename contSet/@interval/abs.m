function res = abs(obj)
% abs - returns the absolute value of an interval
%
% Syntax:  
%    res = abs(obj) 
%
% Inputs:
%    obj - interval object
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
% See also: ---

% Author:       Matthias Althoff
% Written:      26-June-2015
% Last update:  14-February-2015
%               12-October-2015
% Last revision:---

%------------- BEGIN CODE --------------

%init res
res = interval();

%separate computation for scalar and matrix case for efficient computation

%scalar case
if isscalar(obj)
    if obj.sup < 0
        res.inf = abs(obj.sup);
        res.sup = abs(obj.inf);
    elseif obj.inf > 0
        res = obj;
    else
        res.inf = 0;
        res.sup = max(abs(obj.inf), abs(obj.sup));
    end
    
%matrix case
else
    % obj.inf > 0 case
    res = obj;
    
    % find negative indices (if infimum is greater than zero, the absolute value has no effect)
    ind = find(obj.inf < 0 & obj.sup > 0);  % For [-a, +b] case
    ind1 = find(obj.inf < 0 & obj.sup <= 0); % For [-a, -b] case
    
    res.sup(ind) = max(abs(obj.inf(ind)), abs(obj.sup(ind))); % order of computation matters
    res.inf(ind) = abs(0*obj.inf(ind));
    
    res.sup(ind1) = abs(obj.inf(ind1));
    res.inf(ind1) = abs(obj.sup(ind1));
end

%------------- END OF CODE --------------