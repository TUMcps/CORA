function display(loc)
% display - Displays a location object on the command window
%
% Syntax:
%    display(loc)
%
% Inputs:
%    loc - location object
%
% Outputs:
%    ---
%
% Example:
%    inv = interval([-2;-1],[1;2]);
%    trans = transition(conHyperplane([2 0],1),...
%               struct('A',eye(2),'c',zeros(2,1)),2);
%    sys = linearSys([1 -3; -2 1],1);
%    loc = location(inv,trans,sys)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff, Mark Wetzlinger
% Written:       06-November-2007
% Last update:   18-June-2022 (MW, update displayed information)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check if called from display hybridAutomaton
st = dbstack("-completenames");
callFromHybridDisplay = false;
if length(st) >= 2 && contains(st(2).file,['@hybridAutomaton' filesep 'display'])
    callFromHybridDisplay = true;
end

if ~callFromHybridDisplay
    fprintf(newline);
    disp([inputname(1), ' =']);
    fprintf(newline);
end

% array of location objects
if length(loc) > 1
    disp("  " + length(loc) + "x1 location array");
    fprintf(newline);
    return
end

% display name
if strcmp(loc.name,'location')
    disp("  Name: '" + loc.name + "' (default)");
else
    disp("  Name: '" + loc.name + "'");
end

% invariant
if isnumeric(loc.invariant) && isempty(loc.invariant)
    disp("  Invariant: []");
else
    disp("  Invariant: " + class(loc.invariant) + ...
        " (dimension: " + dim(loc.invariant) + ")");
end


% transitions
if isempty(loc.transition)
    disp("  Number of transitions: 0");
else
    transStr = "  Number of transitions: " + length(loc.transition) + " (";
    if isscalar(loc.transition(1).target)
        targetLoc = arrayfun(@(x) x.target,loc.transition,'UniformOutput',true);
        % different grammar...
        if length(loc.transition) == 1
            addString = "target location: ";
        else
            addString = "target locations: ";
        end
        % loop over targets of transition and synchronization labels
        temp = [];
        for i=1:length(loc.transition)
            syncLabel = loc.transition(i).syncLabel;
            if isempty(syncLabel)
                temp = [temp string(targetLoc(i))];
            else
                temp = [temp string(targetLoc(i)) + " ('" + syncLabel + "')"];
            end
        end
        
    else
        % location from location product of parallel hybrid automaton
        temp = [];
        for i=1:length(loc.transition)
            syncLabel = loc.transition(i).syncLabel;
            if isempty(syncLabel)
                temp = [temp "[" + strjoin(string(loc.transition(i).target),",") + "]"];
            else
                temp = [temp "[" + strjoin(string(loc.transition(i).target),",") + "]" ...
                    + " ('" + syncLabel + "')"];
            end
        end
        addString = "target locations: ";
    end
    % extend first entry by additional string
    temp(1) = addString + temp(1);
    % extend last entry by parenthesis
    temp(end) = temp(end) + ")";
    % display transition
    dispUpToLength(temp,100,transStr);
end

% dynamics
disp("  Dynamics: " + class(loc.contDynamics) + ...
    " (state dim.: " + loc.contDynamics.dim + ...
    ", input dim.: " + loc.contDynamics.nrOfInputs + ...
    ", output dim.: " + loc.contDynamics.nrOfOutputs + ")");

if ~callFromHybridDisplay
    fprintf(newline);
end

% ------------------------------ END OF CODE ------------------------------
