function res = test_hybridAutomaton_isequal
% test_hybridAutomaton_isequal - test function for isequal
%
% Syntax:
%    res = test_hybridAutomaton_isequal
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

% compare empty automata
res = isequal(hybridAutomaton(),hybridAutomaton());

% invariant
inv_2D_poly = polytope([-1,0],0);
inv_2D_int = interval([-4;-2],[4;8]);
inv_3D_int = interval([-1;0;0],[1;2;3]);

% transitions
c = [-1;0]; d = 0; C = [0,1]; D = 0;
guard = conHyperplane(c,d,C,D);
reset = struct('A',[1,0;0,-0.75;0,-1],'c',[0;0;0]);
trans_2D_3D = transition(guard,reset,2);

reset = struct('A',[1,0;0,-0.75],'c',[0;0]);
trans_2D_2D = transition(guard,reset,1);

clear reset
reset.f = @(x,u) [x(1);-x(2)*x(1);sqrt(x(1))];
trans_2D_3D_nonlin = transition(guard,reset,2);

c = [-1;0;0]; d = 0; C = [0,0,1]; D = 0;
guard = conHyperplane(c,d,C,D);
reset = struct('A',zeros(2,3),'c',zeros(2,1));
trans_3D_2D = transition(guard,reset,1);

% flow equation
dynamics_2D_lin = linearSys([0,1;0,0],[0;0],[0;-9.81]);
dynamics_2D_nonlin = nonlinearSys(@(x,u) [x(2); -x(1)*x(2)]);
dynamics_3D_lin = linearSys([0,1,0;0,0,1;0,0,1],1);

% define locations
loc_1 = location(inv_2D_poly,trans_2D_3D,dynamics_2D_lin);
loc_2 = location(inv_3D_int,trans_3D_2D,dynamics_3D_lin);
loc_3 = location(inv_2D_poly,trans_2D_2D,dynamics_2D_lin);
loc_4 = location(inv_2D_int,trans_2D_2D,dynamics_2D_lin);
loc_5 = location(inv_2D_int,trans_2D_2D,dynamics_2D_nonlin);
loc_6 = location(inv_2D_int,trans_2D_3D_nonlin,dynamics_2D_nonlin);

% same hybrid automaton
HA1 = hybridAutomaton([loc_1;loc_2]);
res(end+1,1) = isequal(HA1,HA1);

% same array of hybrid automata
temp = isequal([HA1;HA1],[HA1;HA1]);
res(end+1,1) = all(size(temp) == [2,1]);
res(end+1,1) = all(temp);

% different number of locations
HA2 = hybridAutomaton(loc_3);
res(end+1,1) = ~isequal(HA1,HA2);

% different invariant
HA3 = hybridAutomaton(loc_4);
res(end+1,1) = ~isequal(HA2,HA3);

% different flow
HA4 = hybridAutomaton(loc_5);
res(end+1,1) = ~isequal(HA2,HA4);

% different reset function
HA5 = hybridAutomaton([loc_6;loc_2]);
res(end+1,1) = ~isequal(HA2,HA5);

% different length of arrays
res(end+1,1) = ~isequal([HA1;HA2],HA3);

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
