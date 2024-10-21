function res = test_location_isemptyobject
% test_location_isemptyobject - test function for isemptyobject
%
% Syntax:
%    res = test_location_isemptyobject
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       16-May-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty location
assert(isemptyobject(location()));

% non-empty location
inv = polytope([-1,0],0);
guard = polytope([0,1],0,[-1,0],0);
reset = linearReset([1,0;0,-0.75]);
trans = transition(guard,reset,1);
dynamics = linearSys([0,1;0,0],[0;0],[0;-9.81]);

assert(~isemptyobject(location(inv,trans,dynamics)));
assert(all(isemptyobject([location(),location(inv,trans,dynamics)]) == [true false]));

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
