function val = get(obj, propName)
% get - Retrieve object data from obj
%
% Syntax:  
%    val = get(obj, propName)
%
% Properties:
%    guard - guard of the transition
%    target - target location of transition

% Author:       Matthias Althoff
% Written:      03-May-2007
% Last update:  23-Aug-2013
% Last revision:---

%------------- BEGIN CODE --------------

switch propName
    case 'reset'
        val = obj.reset;
    case 'guard'
        val = obj.guard; 
    case 'target'
        val = obj.target;
    case 'equations'
        val = get(obj.guard,'equations');
    otherwise
        error([propName,' is not a valid asset property'])
end

%------------- END OF CODE --------------