function han = plotAsGraph(R)
% plotAsGraph - plot branches of reachable set as graph, branches with only
%    one time-point solution (due to an instant outgoing transition) are
%    marked in red
%
% Syntax:
%    han = plotAsGraph(R)
%
% Inputs:
%    R - reachSet object
%
% Outputs:
%    han - handle to graph object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       07-June-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if length(R) == 1

    % plot graph
    han = plot(graph(1,1,'omitselfloops'));

else

    % convert reachSet object to directed graph
    edgeStart = [R.parent];
    edgeStart = edgeStart(2:end);
    edgeEnd = 2:length(R);
    
    % instantiate directed graph
    G = graph(edgeStart,edgeEnd);

    % node label is location, rewrite as string for plot
    nodeLabel = "[" + string(num2str([R.loc]')) + "]";
    
    % plot graph
    han = plot(G,'Layout','layered','NodeLabel',nodeLabel);

    % find instant transitions (reachSet object with only one time-point
    % solution) and plot these nodes in a different color
    instantTrans = arrayfun(@(x) isempty(x.timeInterval) ...
        && isscalar(x.timePoint.set),R,'UniformOutput',true);

    % use 0 for non-instant transitions and 1 for instant transitions
    han.NodeCData = instantTrans;
    % set color map: blue/red for non-instant/instant transitions
    colormap([0 0 1; 1 0 0]);
    drawnow;

end

% return handle only if desired
if nargout == 0
    clear han;
end

% ------------------------------ END OF CODE ------------------------------
