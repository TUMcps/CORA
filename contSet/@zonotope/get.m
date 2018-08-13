function val = get(obj, propName)
% get - Retrieve object data from obj
%
% Syntax:  
%    val = get(obj, propName)
%
% Properties:
%    Z - zonotope matrix

% Author: Matthias Althoff
% Written: 30-September-2006 
% Last update: 23-March-2007
% Last revision: ---

%------------- BEGIN CODE --------------

switch propName
    case 'Z'
        val = obj.Z;   
    case 'equations'
        val= obj.halfspace.equations;
otherwise
    error([propName,' is not a valid asset property'])
end

%------------- END OF CODE --------------