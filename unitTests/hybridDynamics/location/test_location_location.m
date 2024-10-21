function res = test_location_location
% test_location_location - unit test for constructor of the class location
%
% Syntax:
%    res = test_location_location
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
% Written:       26-November-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty object
loc = location();

% name of the location
name = 'S1';

% invariant
inv = polytope([-1,0],0);

% transition
guard = polytope([0 1],0,[-1 0],0);

% reset function
reset = linearReset([1,0;0,-0.75]);

% transition
trans = transition(guard,reset,2);

% flow equation
dynamics = linearSys([0,1;0,0],[0;0],[0;-9.81]);

% define location
loc = location(name,inv,trans,dynamics);

% check if values have been assigned correctly
assert(strcmp(loc.name,name));
assert(isequal(loc.invariant,inv));
assert(isequal(loc.transition,trans));
assert(isequal(loc.contDynamics,linearSys([0,1;0,0],[0;0],[0;-9.81])));

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
