function res = plus(summand1,summand2)
% plus - Overloaded '+' operator for intervals
%
% Syntax:  
%    res = plus(summand1,summand2)
%
% Inputs:
%    summand1 - interval (for computational efficiency, no single value
%               considered; does not require type checking)
%    summand2 - interval (for computational efficiency, no single value
%               considered; does not require type checking)
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
% Written:      19-June-2015
% Last update:  23-June-2015
%               10-August-2016
%               24-August-2016
%               05-May-2020 (MW, standardized error message)
% Last revision:---

%------------- BEGIN CODE --------------

% determine the interval object
if isa(summand1,'interval')
    res = summand1;
    summand = summand2;
elseif isa(summand2,'interval')
    res = summand2;
    summand = summand1; 
end

% different cases depending on the class of the summand
if isa(summand,'interval')

    res.inf = res.inf + summand.inf;
    res.sup = res.sup + summand.sup;

elseif isnumeric(summand)

    res.inf = res.inf + summand;
    res.sup = res.sup + summand;    

elseif isa(summand,'zonotope') || isa(summand,'conZonotope') || ...
       isa(summand,'zonoBundle') || isa(summand,'polyZonotope') || ...
       isa(summand,'mptPolytope') || isa(summand,'conPolyZono')

    res = summand + res;

else

    % throw error for given arguments
    error(noops(summand1,summand2));
end

%------------- END OF CODE --------------