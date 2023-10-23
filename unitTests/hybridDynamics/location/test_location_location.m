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

res = [];

% empty object
loc = location();

% name of the location
name = 'S1';

% invariant
inv = polytope([-1,0],0);

% transition
c = [-1;0]; d = 0; C = [0,1]; D = 0;
guard = conHyperplane(c,d,C,D);

% reset function
reset = struct('A',[1,0;0,-0.75],'c',[0;0]);

% transition
trans = transition(guard,reset,2);

% flow equation
dynamics = linearSys([0,1;0,0],[0;0],[0;-9.81]);

% define location
loc = location(name,inv,trans,dynamics);

% check if values have been assigned correctly
res(end+1,1) = strcmp(loc.name,name);
res(end+1,1) = isequal(loc.invariant,inv);
res(end+1,1) = isequal(loc.transition,trans);
res(end+1,1) = isequal(loc.contDynamics,linearSys([0,1;0,0],[0;0],[0;-9.81]));

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
