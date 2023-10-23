function PHA = roomHeatingParallel()
% roomHeatingParallel - room heating benchmark described in Sec. 2.3 in [1]
%                       with two rooms represented as a parallel hybrid
%                       automaton
%
% Syntax:  
%    PHA = roomHeatingParallel()
%
% Inputs:
%    ---
%
% Outputs:
%    PHA - parallelHybridAutomaton object
%
% References:
%   [1] A. Fehnker and F. Ivancic. "Benchmarks for Hybrid Systems 
%       Verification", HSCC 2004

% Author:       Niklas Kochdumper
% Written:      26-June-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------


% Test parallelHybridAutomaton

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

linSys = linearSys(A,B,c,C);

inv = polytope(1,T_off);

guard = conHyperplane(1,T_off);

reset.A = 1;
reset.c = 0;

trans = transition(guard,reset,2);

loc = location('on',inv,trans,linSys);

% Location 2 : Heating off
A = -(a1 + b1);
B = [b1,a1];
C = 1;

linSys = linearSys(A,B,[],C);

inv = polytope(-1,-T_on);

guard = conHyperplane(1,T_on);

reset.A = 1;
reset.c = 0;

trans = transition(guard,reset,1);

loc(2) = location('off',inv,trans,linSys);

% Hybrid automaton
HA1 = hybridAutomaton(loc);



% Room 2 - HA 2 -----------------------------------------------------------

% Location 1 : Heating on
A = -(a2 + b2);
B = [b2,a2];
c = c2;
C = 1;

linSys = linearSys(A,B,c,C);

inv = polytope(1,T_off);

guard = conHyperplane(1,T_off);

reset.A = 1;
reset.c = 0;

trans = transition(guard,reset,2);

loc(1) = location('on',inv,trans,linSys);

% Location 2 : Heating off
A = -(a2 + b2);
B = [b2,a2];
C = 1;

linSys = linearSys(A,B,[],C);

inv = polytope(-1,-T_on);

guard = conHyperplane(1,T_on);

reset.A = 1;
reset.c = 0;

trans = transition(guard,reset,1);

loc(2) = location('off',inv,trans,linSys);

% Hybrid automaton
HA2 = hybridAutomaton(loc);


% Parallel Hybrid Automaton -----------------------------------------------

% components
comp = [HA1;HA2];

% connections between the components
inputBinds{1} = [0 1; ...   % first global input
                 2 1];      % first output of component 2
             
inputBinds{2} = [0 1; ...   % first global input
                 1 1];      % first output of component 1
   
% parallel hybrid automaton
PHA = parallelHybridAutomaton(comp,inputBinds);

end

%------------- END OF CODE --------------