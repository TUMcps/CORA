function viewGraph(HA)
% viewGraph - Plots the discrete structure of a hybrid automaton as a graph
%
% Syntax:
%    viewGraph(HA)
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
%    guard = polytope([0,1],0,[-1,0],0);
%    reset = linearReset([1,0;0,-0.75]);
%    trans(1) = transition(guard,reset,1);
% 
%    % flow equation
%    dynamics = linearSys([0,1;0,0],[0;0],[0;-9.81]);
% 
%    % define location
%    loc(1) = location('S1',inv,trans,dynamics);
% 
%    % instantiate hybrid automaton (and display)
%    HA = hybridAutomaton('HA',loc);
%
%    % plot graph for the hybrid automaton
%    viewGraph(HA);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: display

% Authors:       Niklas Kochdumper
% Written:       15-September-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % loop over all transitions
    s = []; t = [];

    for i = 1:length(HA.location)
        for j = 1:length(HA.location(i).transition)
            s = [s,i];
            t = [t,HA.location(i).transition(j).target];
        end
    end

    G = digraph(s,t);

    % catch special case where the final locations do not have transitions
    m = max([s,t]);

    if m < length(HA.location)
        G = G.addnode(length(HA.location)-m);
    end

    % extract names of the locations
    nodeLabels = cell(length(HA.location),1);

    for i = 1:length(HA.location)
        if ~strcmp(HA.location(i).name,'location')
            nodeLabels{i} = HA.location(i).name;
        else
            nodeLabels{i} = num2str(i);
        end
    end

    % visualize the graph
    plot(G,'NodeLabel',nodeLabels);
end

% ------------------------------ END OF CODE ------------------------------
