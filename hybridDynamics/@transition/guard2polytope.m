function obj = guard2polytope(obj)
% guard2polytope - convert the guard set to a mptPolytope
%
% Syntax:  
%    obj = guard2polytope(obj)
%
% Inputs:
%    obj - transition object
%
% Outputs:
%    obj - modified transition object
%
% Other m-files required: not specified
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Niklas Kochdumper
% Written:      16-May-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
    
    if ~isa(obj.guard,'levelSet') && ~isa(obj.guard,'mptPolytope')
        obj.guard = mptPolytope(obj.guard);
    end
end

%------------- END OF CODE --------------