function pHA = roomHeatingParallel()
% roomHeatingParallel - room heating benchmark described in Sec. 2.3 in [1]
%    with two rooms represented as a parallel hybrid automaton
%
% Syntax:
%    pHA = roomHeatingParallel()
%
% Inputs:
%    -
%
% Outputs:
%    pHA - parallelHybridAutomaton object
%
% References:
%   [1] A. Fehnker and F. Ivancic. "Benchmarks for Hybrid Systems 
%       Verification", HSCC 2004
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Niklas Kochdumper
% Written:       26-June-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parameter room 1
a1 = 0.5; 
b1 = 0.4; 
c1 = 6;

% parameter room 2
a2 = 0.5;
b2 = 0.3;
c2 = 7;

% temperatures for switching heaters on and off
T_off = 21;
T_on = 20;


% Room 1 - HA 1 -----------------------------------------------------------

% Location 1 : Heating on
A = -(a1 + b1);
B = [b1,a1];
c = c1;
C = 1;

linsys = linearSys(A,B,c,C);

inv = polytope(1,T_off);

guard = polytope([],[],1,T_off);

reset = linearReset(1,0,0);

trans = transition(guard,reset,2);

loc = location('on',inv,trans,linsys);

% Location 2 : Heating off
A = -(a1 + b1);
B = [b1,a1];
C = 1;

linsys = linearSys(A,B,[],C);

inv = polytope(-1,-T_on);

guard = polytope([],[],1,T_on);

reset = linearReset(1,0,0);

trans = transition(guard,reset,1);

loc(2) = location('off',inv,trans,linsys);

% Hybrid automaton
HA1 = hybridAutomaton('roomHeatingParallel1',loc);


% Room 2 - HA 2 -----------------------------------------------------------

% Location 1 : Heating on
A = -(a2 + b2);
B = [b2,a2];
c = c2;
C = 1;

linsys = linearSys(A,B,c,C);

inv = polytope(1,T_off);

guard = polytope([],[],1,T_off);

reset = linearReset(1,0,0);

trans = transition(guard,reset,2);

loc(1) = location('on',inv,trans,linsys);

% Location 2 : Heating off
A = -(a2 + b2);
B = [b2,a2];
C = 1;

linsys = linearSys(A,B,[],C);

inv = polytope(-1,-T_on);

guard = polytope([],[],1,T_on);

reset = linearReset(1,0,0);

trans = transition(guard,reset,1);

loc(2) = location('off',inv,trans,linsys);

% Hybrid automaton
HA2 = hybridAutomaton('roomHeatingParallel2',loc);


% Parallel Hybrid Automaton -----------------------------------------------

% components
comp = [HA1;HA2];

% connections between the components
inputBinds{1} = [0 1; ...   % first global input
                 2 1];      % first output of component 2
             
inputBinds{2} = [0 1; ...   % first global input
                 1 1];      % first output of component 1
   
% parallel hybrid automaton
pHA = parallelHybridAutomaton('roomHeatingParallel',comp,inputBinds);

% ------------------------------ END OF CODE ------------------------------
