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
%    inv = polytope([-1,0],0);
% 
%    % transition
%    c = [-1;0]; d = 0; C = [0,1]; D = 0;
%    guard = conHyperplane(c,d,C,D);
%    reset = struct('A',[1,0;0,-0.75],'c',[0;0]);
%    trans(1) = transition(guard,reset,1);
% 
%    % flow equation
%    dynamics = linearSys([0,1;0,0],[0;0],[0;-9.81]);
% 
%    % define location
%    loc(1) = location('S1',inv,trans,dynamics);
% 
%    % instantiate hybrid automaton (and display)
%    HA = hybridAutomaton(loc)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       18-June-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

fprintf(newline);

disp([inputname(1), ' =']);

fprintf(newline);

if isemptyobject(HA)

    dispEmptyObj(HA,inputname(1));

elseif length(HA) > 1
    
    disp("  " + length(HA) + "x1 hybridAutomaton array");

else

    % number of locations
    numLoc = length(HA.location);
    
    % loop over locations
    for i=1:numLoc
        
        % number of location
        disp("Location " + i + " of " + numLoc + ":");
        display(HA.location(i));
    
    end
end

fprintf(newline);

% ------------------------------ END OF CODE ------------------------------
