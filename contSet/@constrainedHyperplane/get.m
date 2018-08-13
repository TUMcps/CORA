function val = get(obj, propName)
% get - Retrieve object data from obj
%
% Syntax:  
%    val = get(obj, propName)
%
% Properties:
%    equations - number of halfspace equations

% Author:       Matthias Althoff
% Written:      10-August-2011
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

switch propName 
    case 'equations'
        val = 1; %return number of halfspace equations
    otherwise
        error([propName,' is not a valid asset property'])
end

%------------- END OF CODE --------------