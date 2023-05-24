function trans = guard2polytope(trans)
% guard2polytope - convert the guard set of a transition to a mptPolytope
%
% Syntax:  
%    trans = guard2polytope(trans)
%
% Inputs:
%    trans - transition object
%
% Outputs:
%    trans - modified transition object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Niklas Kochdumper
% Written:      16-May-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

if ~isa(trans.guard,'levelSet') && ~isa(trans.guard,'mptPolytope')
    trans.guard = mptPolytope(trans.guard);
end

%------------- END OF CODE --------------