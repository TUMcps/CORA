function display(pHA)
% display - Displays the properties of a parallelHybridAutomaton object on
%    the command window
%
% Syntax:
%    display(pHA)
%
% Inputs:
%    pHA - hybridAutomaton object
%
% Outputs:
%    -
%
% Example: 
%    % parameters
%    a1 = 0.5; b1 = 0.4; c1 = 6;
%    a2 = 0.5; b2 = 0.3; c2 = 7;
%    T_off = 21; T_on = 20;
% 
%    % first component, first location
%    inv = polytope(1,T_off);
%    guard = conHyperplane(1,T_off);
%    reset = struct('A',1,'c',0);
%    trans = transition(guard,reset,2);
%    linSys = linearSys(-(a1 + b1),[b1,a1],c1,1);
%    loc(1) = location('on',inv,trans,linSys);
% 
%    % first component, second location
%    inv = polytope(-1,-T_on);
%    guard = conHyperplane(1,T_on);
%    reset = struct('A',1,'c',0);
%    trans = transition(guard,reset,1);
%    linSys = linearSys(-(a1 + b1),[b1,a1],[],1);
%    loc(2) = location('off',inv,trans,linSys);
% 
%    HA1 = hybridAutomaton(loc);
% 
%    % second component, first location
%    inv = polytope(1,T_off);
%    guard = conHyperplane(1,T_off);
%    reset = struct('A',1,'c',0);
%    trans = transition(guard,reset,2);
%    linSys = linearSys(-(a2 + b2),[b2,a2],c2,1);
%    loc(1) = location('on',inv,trans,linSys);
% 
%    % second component, second location
%    inv = polytope(-1,-T_on);
%    guard = conHyperplane(1,T_on);
%    reset = struct('A',1,'c',0);
%    trans = transition(guard,reset,1);
%    linSys = linearSys(-(a2 + b2),[b2,a2],[],1);
%    loc(2) = location('off',inv,trans,linSys);
% 
%    HA2 = hybridAutomaton(loc);
% 
%    % parallel hybrid automaton
%    components = [HA1;HA2];
%    inputBinds{1} = [0 1; ...   % first global input
%                     2 1];      % first output of component 2
%    inputBinds{2} = [0 1; ...   % first global input
%                     1 1];      % first output of component 1
%    pHA = parallelHybridAutomaton(components,inputBinds);
%    display(pHA);
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

if all(isemptyobject(pHA))

    dispEmptyObj(pHA,inputname(1));
    return

elseif length(pHA) > 1

    disp("  " + length(pHA) + "x1 parallelHybridAutomaton object");
    fprintf(newline);
    return

end

% number of components
numComp = length(pHA.components);

% loop over components
for i=1:numComp
    
    % number of component
    disp("Component " + i + " of "+ numComp + ":");

    % number of locations
    numLoc = length(pHA.components(i).location);
    locString = "   Number of locations: " + numLoc + " (";
    % gather names of locations
    if all(arrayfun(@(x) strcmp(x.name,'location'),...
            pHA.components(i).location,'UniformOutput',true))
        % no location has an actual name (all default names)
        temp = "no names";
    else
        temp = [];
        for j=1:numLoc
            nameLoc = pHA.components(i).location(j).name;
            if strcmp(nameLoc,'location')
                % default name
                temp = [temp "(no name)"];
            else
                temp = [temp "'" + string(nameLoc) + "'"];
            end
        end
    end
    % extend last entry by closing parenthesis
    temp(end) = temp(end) + ")";
    dispUpToLength(temp,100,locString);

    
    % transitions in location
    transString = [];
    for j=1:numLoc
        numTrans = length(pHA.components(i).location(j).transition);
        target = []; syncLabel = {};
        % read out target locations and synchronization labels
        for k=1:numTrans
            target(end+1,1) = pHA.components(i).location(j).transition(k).target;
            syncLabel{end+1,1} = pHA.components(i).location(j).transition(k).syncLabel;
        end
        % go over all transitions
        while ~isempty(target)
            % number of transition with same target location
            numTarget = nnz(target == target(1));
            if numTarget == 1
                if ~isempty(syncLabel{1})
                    transString = [transString "loc" + j + " -> loc" + target(1) ...
                                    + " ('" + syncLabel{1} + "')"];
                else
                    transString = [transString "loc" + j + " -> loc" + target(1)];
                end
            else % numTarget > 1
                temp = syncLabel(target == target(1));
                temp = temp(cellfun(@(x)~isempty(x),temp,'UniformOutput',true));
                if ~isempty(temp)
                    transString = [transString "loc" + j + " -> loc" + target(1) ...
                        + " (" + numTarget + " times, incl. labels: '" + strjoin(temp,"','") + "')"];
                else
                    transString = [transString "loc" + j + " -> loc" + target(1) ...
                        + " (" + numTarget + "times)"];
                end
            end
            % remove all transitions with same target location
            syncLabel = syncLabel(target ~= target(1));
            target = target(target ~= target(1));
        end
    end
    % check if string needs to be split in separate lines
    dispUpToLength(transString,100,"   Transitions: ");

    % state dimension, input dimension (including input binds), output
    % dimension (has to be equal over all locations of a component)
    disp("   State dimension: " + pHA.components(i).location(1).contDynamics.dim);

    % input binds
    if isempty(pHA.bindsInputs{i})
        disp("   Input dimension: 0");
    else
        inpString = "   Input dimension: " + ...
            pHA.components(i).location(1).contDynamics.nrOfInputs + " (";
        temp = [];
        % loop over each input
        for j=1:pHA.components(i).location(1).contDynamics.nrOfInputs
            origin_comp = pHA.bindsInputs{i}(j,1);
            origin_idx = pHA.bindsInputs{i}(j,2);
            if origin_comp == 0
                % global
                temp = [temp "u" + j + " = u" + origin_idx + "(global)"];
            else
                % output from other component
                temp = [temp "u" + j + " <- y" + origin_idx + " of comp.#" + origin_comp];
            end
        end
        % extend last entry by parenthesis
        temp(end) = temp(end) + ")";
        dispUpToLength(temp,100,inpString);
    end
    
    % output dimension
    disp("   Output dimension: " + pHA.components(i).location(1).contDynamics.nrOfOutputs);

    fprintf(newline);

end

% ------------------------------ END OF CODE ------------------------------
