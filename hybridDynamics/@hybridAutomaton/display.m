function display(HA)
% display - Displays the properties of a hybridAutomaton object on the
%    command window
%
% Syntax:  
%    display(HA)
%
% Inputs:
%    HA - hybridAutomaton object
%
% Outputs:
%    -
%
% Example: 
%    % invariant
%    polyOpt = struct('A',[-1,0],'b',0);
%    inv = mptPolytope(polyOpt);
% 
%    % transition
%    c = [-1;0]; d = 0; C = [0,1]; D = 0;
%    guard = conHyperplane(c,d,C,D);
%    reset = struct('A',[1,0;0,-0.75],'c',[0;0]);
%    trans{1} = transition(guard,reset,1);
% 
%    % flow equation
%    dynamics = linearSys([0,1;0,0],[0;0],[0;-9.81]);
% 
%    % define location
%    loc{1} = location('S1',inv,trans,dynamics);
% 
%    % instantiate hybrid automaton (and display)
%    HA = hybridAutomaton(loc)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Mark Wetzlinger
% Written:      18-June-2022
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

fprintf(newline);

disp([inputname(1), ' =']);

fprintf(newline);

% number of locations
numLoc = length(HA.location);

% loop over locations
for i=1:numLoc
    
    % number of location
    disp("Location " + i + " of " + numLoc + ":");
    display(HA.location{i});

end

fprintf(newline);

%------------- END OF CODE --------------