function val = get(obj, propName)
% get - Retrieve object data from obj
%
% Syntax:  
%    val = get(obj, propName)
%
% Properties:
%    taylor - taylor series related properties

% Author: Matthias Althoff
% Written: 30-October-2007
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------

switch propName
    case 'taylor'
        val = obj.taylor;   
    case 'A'
        val = obj.A; 
otherwise
    error([propName,' is not a valid asset property'])
end

%------------- END OF CODE --------------