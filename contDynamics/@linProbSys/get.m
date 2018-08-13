function val = get(obj, propName)
% get - Retrieve object data from obj
%
% Syntax:  
%    val = get(obj, propName)
%
% Properties:
%    Z - zonotope matrix

% Author:       Matthias Althoff
% Written:      28-September-2007
% Last update:  03-September-2009
% Last revision:---

%------------- BEGIN CODE --------------

switch propName
    case 'taylor'
        val = obj.taylor;    
    case 'eAt'
        val = obj.taylor.eAt;   
otherwise
    error([propName,' is not a valid asset property'])
end

%------------- END OF CODE --------------