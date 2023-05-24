function res = test_hybridAutomaton_dim
% test_hybridAutomaton_dim - test function for dim
%
% Syntax:  
%    res = test_hybridAutomaton_dim
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

% Author:       Mark Wetzlinger
% Written:      23-May-2023
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% generate simple automaton
inv1 = mptPolytope([-1 1]/sqrt(2),0);
inv2 = mptPolytope([1 -1 0]/sqrt(2),0);

guard1 = conHyperplane([-1 1]/sqrt(2),0);
guard2 = conHyperplane([-1 1 0]/sqrt(2),0);
reset1 = struct('A',[1 0; 0 1; 0 0],'c',[-1;1;1]);
trans1 = transition(guard1,reset1,2);

reset2 = struct('A',[1 0 0; 0 1 0],'c',[1;-1]);
trans2 = transition(guard2,reset2,1);

dyn1 = linearSys([0 1; -1 0],1);
dyn2 = linearSys([0 -1 0; 1 0 0; 0 0 0],1);

loc1 = location('clockwise',inv1,trans1,dyn1);
loc2 = location('counter-clockwise',inv2,trans2,dyn2);
HA = hybridAutomaton([loc1;loc2]);

% dimension
n = dim(HA);

res = length(n) == 2 && isnumeric(n);
res(end+1,1) = n(1) == 2 && n(2) == 3;

% array of hybrid automata
n = dim([HA;HA]);

res(end+1,1) = iscell(n) && length(n) == 2;
res(end+1,1) = n{1}(1) == 2 && n{1}(2) == 3 ...
    && n{2}(1) == 2 && n{2}(2) == 3;


% combine results
res = all(res);

%------------- END OF CODE --------------
