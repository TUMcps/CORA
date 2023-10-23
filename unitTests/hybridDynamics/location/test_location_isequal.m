function res = test_location_isequal
% test_location_isequal - test function for isequal
%
% Syntax:
%    res = test_location_isequal
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

% empty locations
res = isequal(location(),location());

% invariant
inv1 = polytope([-1,0],0);
inv2 = polytope([-1,0],0.5);

% transition
guard = conHyperplane([-1;0],0,[0,1],0);

% reset function
reset1 = struct('A',[1,0;0,-0.75],'c',[0;0]);
reset2 = struct('A',[1,0;0,-0.75],'c',[1;0]);

% transition
trans1 = transition(guard,reset1,1);
trans2 = transition(guard,reset1,2);
trans3 = transition(guard,reset2,2);

% flow equation
dynamics1 = linearSys([0,1;0,0],[0;0],[0;-9.81]);
dynamics2 = linearSys([0,1;0,0],[0;0],[0;-9.80]);


% same location
res(end+1,1) = isequal(location(inv1,trans1,dynamics1),...
    location(inv1,trans1,dynamics1));
res(end+1,1) = isequal(location(inv1,[trans1;trans2],dynamics1),...
    location(inv1,[trans2;trans1],dynamics1));

% same array of locations
temp = isequal([location(inv1,trans1,dynamics1),location(inv2,trans2,dynamics2)],...
    [location(inv1,trans1,dynamics1),location(inv2,trans2,dynamics2)]);
res(end+1,1) = all(size(temp) == [1,2]);
res(end+1,1) = all(temp);

% different invariant
res(end+1,1) = ~isequal(location(inv1,trans1,dynamics1),...
    location(inv2,trans1,dynamics1));

% different dynamics
res(end+1,1) = ~isequal(location(inv1,trans1,dynamics1),...
    location(inv1,trans1,dynamics2));

% different transition
res(end+1,1) = ~isequal(location(inv1,trans1,dynamics1),...
    location(inv1,trans2,dynamics1));

% different transition
res(end+1,1) = ~isequal(location(inv1,trans1,dynamics1),...
    location(inv1,[trans1;trans3],dynamics1));

% different array length
res(end+1,1) = ~isequal(location(inv1,trans1,dynamics1),...
    [location(inv1,trans1,dynamics1),location(inv2,trans2,dynamics2)]);

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
