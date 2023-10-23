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
res = isemptyobject(location());

% non-empty location
inv = polytope([-1,0],0);
guard = conHyperplane([-1;0],0,[0,1],0);
reset = struct('A',[1,0;0,-0.75],'c',[0;0]);
trans = transition(guard,reset,1);
dynamics = linearSys([0,1;0,0],[0;0],[0;-9.81]);

res(end+1,1) = ~isemptyobject(location(inv,trans,dynamics));
res(end+1,1) = all(isemptyobject([location(),location(inv,trans,dynamics)]) == [true false]);

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
