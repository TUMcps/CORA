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
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Mark Wetzlinger
% Written:      26-November-2022
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% assume true
res = true;

% name of the location
name = 'S1';

% invariant
polyOpt = struct('A',[-1,0],'b',0);
inv = mptPolytope(polyOpt);

% transition
c = [-1;0]; d = 0; C = [0,1]; D = 0;
guard = conHyperplane(c,d,C,D);

% reset function
reset = struct('A',[1,0;0,-0.75],'c',[0;0]);

% transition
trans{1} = transition(guard,reset,2);

% flow equation
dynamics = linearSys([0,1;0,0],[0;0],[0;-9.81]);

% define location
loc = location(name,inv,trans,dynamics);

% check if values have been assigned correctly
if ~strcmp(loc.name,name)
    res = false;
end
if ~isequal(loc.invariant,inv)
    res = false;
end
if ~isequal(loc.transition{1},trans{1})
    res = false;
end
if ~isequal(loc.contDynamics,linearSys([0,1;0,0],[0;0],[0;-9.81]))
    res = false;
end

%------------- END OF CODE --------------