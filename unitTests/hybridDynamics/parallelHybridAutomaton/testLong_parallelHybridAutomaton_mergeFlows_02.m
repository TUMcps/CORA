function res = testLong_parallelHybridAutomaton_mergeFlows_02
% testLong_parallelHybridAutomaton_mergeFlows_02 - test function
%    for merging flow equations in the location product; mixed linear/
%    nonlinear flow
%
% Syntax:
%    res = testLong_parallelHybridAutomaton_mergeFlows_02
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
% Written:       22-January-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume true
res = true;

% comp1: x1,x2,     u1,u2,  y1
% comp2: x3,x4,x5,  u3,     y2,y3
syms x1 x2 x3;

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
dynamics{1,2} = linearSys([0 -1; -1 0],[0 1; 1 0],[],[0.05 -0.05]);
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
dynamics{2,1} = linearSys([0 1 -1; 1 0 0; 0 1 0],[0;-1;0],[],[0 0 0.05; 0.05 0.05 0]);
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


% merge flows of location pair [1,1]
mergedFlow = mergeFlows(pHA,{pHA.components(1).location(1).contDynamics,...
    pHA.components(2).location(1).contDynamics},[1;1]);

% instantiate true solution
% comp#1/loc#1: f = @(x,u) [x(2)-u(2); -x(1)^2+u(1)];
%               u(1) = y(1) of comp #2, u(2) = y(2) of comp #2
%               g = @(x,u) x(2)-x(1);
% comp#2/loc#1: f = [0 1 -1; 1 0 0; 0 1 0]*x + [0;-1;0]*u
%               u(1) = y(1) of comp #1
%               g = [0 0 0.05; 0.05 0.05 0]*x;
f = @(x,u) [x(2)-0.05*x(3)-0.05*x(4); ...
           -x(1)^2+0.05*x(5); ...
            x(4)-x(5); ...
            x(3)-x(2)+x(1); ...
            x(4)];
mergedFlow_ = nonlinearSys('location11',f);

% compare solutions
res = isequal(mergedFlow,mergedFlow_,1e-14);

% merge flows of location pair [2,2]
mergedFlow = mergeFlows(pHA,{pHA.components(1).location(2).contDynamics,...
    pHA.components(2).location(2).contDynamics},[2;2]);

% instantiate true solution
% comp#1/loc#2: f = [0 -1; -1 0]*x + [0 1; 1 0]*u
%               u(1) = y(1) of comp #2, u(2) = y(2) of comp #2
%               g = [0.05 -0.05]*x;
% comp#2/loc#2: f = [-x(2); x(1)^2-u(1); x(2)*x(3)];
%               u(1) = y(1) of comp #1
%               g = @(x,u) [x(1)-x(2); x(3)];
f = @(x,u) [-x(2)+x(5); ...
            -x(1)+x(3)-x(4); ...
            -x(4); ...
             x(3)^2-0.05*x(1)+0.05*x(2); ...
             x(4)*x(5)];
mergedFlow_ = nonlinearSys('location22',f);

% compare solutions
res(end+1,1) = isequal(mergedFlow,mergedFlow_,1e-14);

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
