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
% Written:      14-February-2011
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

switch propName
    case 'name'
        val = obj.name;
    case 'invariant'
        val = obj.invariant; 
    case 'transition'
        val = obj.transition;
    case 'contDynamics'
        val = obj.contDynamics;
    otherwise
        error([propName,' is not a valid asset property'])
end

%------------- END OF CODE --------------