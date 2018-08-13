function val = get(obj, propName)
% get - Retrieve object data from obj
%
% Syntax:  
%    val = get(obj, propName)
%
% Properties:
%    equations - number of halfspace equations

% Author:       Niklas Kochdumper
% Written:      16-May-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

switch propName 
    case 'equations'
        % return number of inequatlity constraints describing the interval
        val = 2*length(obj); 
    otherwise
        error([propName,' is not a valid asset property'])
end

%------------- END OF CODE --------------