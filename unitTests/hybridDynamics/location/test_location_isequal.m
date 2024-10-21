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
assert(isequal(location(),location()));

% invariant
inv1 = polytope([-1,0],0);
inv2 = polytope([-1,0],0.5);

% transition
guard = polytope([0,1],0,[-1,0],0);

% reset function
reset1 = linearReset([1,0;0,-0.75]);
reset2 = linearReset([1,0;0,-0.75],[],[1;0]);

% transition
trans1 = transition(guard,reset1,1);
trans2 = transition(guard,reset1,2);
trans3 = transition(guard,reset2,2);

% flow equation
dynamics1 = linearSys([0,1;0,0],[0;0],[0;-9.81]);
dynamics2 = linearSys([0,1;0,0],[0;0],[0;-9.80]);


% same location
assert(isequal(location(inv1,trans1,dynamics1),location(inv1,trans1,dynamics1)));
assert(isequal(location(inv1,[trans1;trans2],dynamics1),location(inv1,[trans2;trans1],dynamics1)));

% same array of locations
log12 = isequal([location(inv1,trans1,dynamics1),location(inv2,trans2,dynamics2)],...
    [location(inv1,trans1,dynamics1),location(inv2,trans2,dynamics2)]);
assert(all(size(log12) == [1,2]));
assert(all(log12));

% different invariant
assert(~isequal(location(inv1,trans1,dynamics1),...
    location(inv2,trans1,dynamics1)));

% different dynamics
assert(~isequal(location(inv1,trans1,dynamics1),location(inv1,trans1,dynamics2)));

% different transition
assert(~isequal(location(inv1,trans1,dynamics1),location(inv1,trans2,dynamics1)));

% different transition
assert(~isequal(location(inv1,trans1,dynamics1),location(inv1,[trans1;trans3],dynamics1)));

% different array length
assert(~isequal(location(inv1,trans1,dynamics1),...
    [location(inv1,trans1,dynamics1),location(inv2,trans2,dynamics2)]));

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
