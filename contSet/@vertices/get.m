function val = get(obj, propName)
% get - Retrieve object data from obj
%
% Syntax:  
%    val = get(obj, propName)
%
% Properties:
%    property1 - vertices

% Author: Matthias Althoff
% Written: 30-September-2006 
% Last update: 23-March-2007
% Last revision: ---

%------------- BEGIN CODE --------------

switch propName
    case 'V'
        val = obj.V;   
    otherwise
        error([propName,' is not a valid asset property'])
end

%------------- END OF CODE --------------