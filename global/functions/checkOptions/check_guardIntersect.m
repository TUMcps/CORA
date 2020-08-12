function check_guardIntersect(options, obj)
% check_guardIntersect - checks if options.guardIntersect
%  1) exists
%  2) takes an allowed value
%
% Syntax:
%    check_guardIntersect(options, obj)
%
% Inputs:
%    options - options for object
%    obj     - system object
%
% Outputs:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none
%
% References: 
%   -

% Author:       Mark Wetzlinger
% Written:      05-Mar-2019
% Last update:  03-May-2020 (rewriting of error msgs using class(obj))
% Last revision:---

%------------- BEGIN CODE --------------

option = 'guardIntersect';
strct = 'options';
% guardIntersect:
% polytope, conZonotope, conZonotopeFast, zonoGirard or pancake
validGuardIntersect = {'polytope', 'conZonotope','levelSet',...
    'zonoGirard', 'pancake', 'hyperplaneMap', 'nondetGuard'};
if ~isfield(options,option)
    error(printOptionMissing(obj,option,strct));
else
    found = 0;
    for i=1:length(validGuardIntersect)
        if strcmp(options.guardIntersect, validGuardIntersect{i})
            found = i; break;
        end
    end
    if ~found
        error(printOptionOutOfRange(obj,option,strct));
    end
end

end

%------------- END OF CODE --------------

