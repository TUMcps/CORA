function res = testLong_parallelHybridAutomaton_mergeInvariants_01
% testLong_parallelHybridAutomaton_mergeInvariants_01 - test
%    function for merging invariant sets (level sets) in the location product
%
% Syntax:
%    res = testLong_parallelHybridAutomaton_mergeInvariants_01
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
% Written:       15-January-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% comp1: x1,x2,     u1,u2,  y1
% comp2: x3,x4,x5,  u3,     y2,y3
syms x1 x2 x3 x4 x5;

% first component, first location
f = @(x,u) [x(2)-u(2); -x(1)^2+u(1)];
g = @(x,u) x(2)-x(1);
dynamics{1,1} = nonlinearSys(f,g);
vars = sym('x',[2,1]);
eq = -x1 + x2^2 + 3;
compOp = '<=';
inv{1,1} = levelSet(eq,vars,compOp);
compOp = '==';
guard = levelSet(eq,vars,compOp);
reset = struct('f',@(x,u) [x(1)-u(1);x(2)-u(2)]);
trans{1,1} = transition(guard,reset,2);
loc = location('loc1',inv{1,1},trans{1,1},dynamics{1,1});

% first component, second location
f = @(x,u) [-x(2)^2+u(2); x(1)+u(1)];
g = @(x,u) x(1)-x(2);
dynamics{1,2} = nonlinearSys(f,g);
eq = x1 - x2^2 - 3;
compOp = '<=';
inv{1,2} = levelSet(eq,vars,compOp);
compOp = '==';
guard = levelSet(eq,vars,compOp);
reset = struct('f',@(x,u) [x(1)-u(1);-x(2)+u(2)]);
trans{1,2} = transition(guard,reset,1);
loc(2) = location('loc2',inv{1,2},trans{1,2},dynamics{1,2});

% init first hybrid automaton
HA1 = hybridAutomaton(loc);

% second component, first location
f = @(x,u) [x(2); x(1)^2-u(1); -x(2)*x(3)];
g = @(x,u) [x(2)-x(1); x(3)];
dynamics{2,1} = nonlinearSys(f,g);
vars = sym('x',[3,1]);
eq = x1^2 - x2 + x3 - 5;
compOp = '<=';
inv{2,1} = levelSet(eq,vars,compOp);
compOp = '==';
guard = levelSet(eq,vars,compOp);
reset = struct('f',@(x,u) [-x(1)+u(1);x(2);x(3)]);
trans{2,1} = transition(guard,reset,2);
loc(1) = location('loc1',inv{2,1},trans{2,1},dynamics{2,1});

% second component, second location
f = @(x,u) [-x(2); x(1)^2-u(1); x(2)*x(3)];
g = @(x,u) [x(1)-x(2); x(3)];
dynamics{2,2} = nonlinearSys(f,g);
eq = -x1^2 + x2 - x3 + 5;
compOp = '<=';
inv{2,2} = levelSet(eq,vars,compOp);
compOp = '==';
guard = levelSet(eq,vars,compOp);
reset = struct('f',@(x,u) [-x(1)-u(1);x(2);-x(3)]);
trans{2,2} = transition(guard,reset,1);
loc(2) = location('loc2',inv{2,2},trans{2,2},dynamics{2,2});

% init second hybrid automaton
HA2 = hybridAutomaton(loc);

% components and input binds
components = [HA1;HA2];
inputBinds{1} = [2 1; 2 2];
inputBinds{2} = [1 1];

% instantiate parallel hybrid automaton
pHA = parallelHybridAutomaton(components,inputBinds);


% merge invariants of location pair [1,2]

% merge invariant sets
mergedInv = mergeInvariants(pHA,...
    {pHA.components(1).location(1).invariant,...
    pHA.components(2).location(2).invariant},...
    {transition();transition()});

% instantiate true solution
eq1 = -x1 + x2^2 + 3;
eq2 = -x3^2 + x4 - x5 + 5;
vars = sym('x',[5,1]);
mergedInv_ = levelSet([eq1;eq2],vars,{'<=';'<='});

% compare solutions
res = isequal(mergedInv,mergedInv_,1e-14);


% merge invariants of location pair [2,1]

% merge invariant sets
mergedInv = mergeInvariants(pHA,...
    {pHA.components(1).location(2).invariant,...
    pHA.components(2).location(1).invariant}, ...
    {transition();transition()});

% instantiate true solution
eq1 = x1 - x2^2 - 3;
eq2 = x3^2 - x4 + x5 - 5;
vars = sym('x',[5,1]);
mergedInv_ = levelSet([eq1;eq2],vars,{'<=';'<='});

% compare solutions
res(end+1,1) = isequal(mergedInv,mergedInv_,1e-14);

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
