function res = test_parallelHybridAutomaton_display
% test_parallelHybridAutomaton_display - test function for display
%
% Syntax:
%    res = test_parallelHybridAutomaton_display
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
% Written:       08-May-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume true
res = true;

% empty parallel hybrid automaton
parallelHybridAutomaton()

% call model function
roomHeatingParallel()


% init model directly
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
loc(1) = location('loc1',inv{1,1},trans{1,1},dynamics{1,1});

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
pHA = parallelHybridAutomaton(components,inputBinds)


% ------------------------------ END OF CODE ------------------------------
