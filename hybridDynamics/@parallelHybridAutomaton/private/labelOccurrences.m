function allLabels = labelOccurrences(pHA)
% labelOccurrences - create list of label occurences to check whether all
%    labeled transitions are enabled at the same time
%
% Syntax:
%    allLabels = labelOccurrences(pHA)
%
% Inputs:
%    pHA - parallelHybridAutomaton object
%
% Outputs:
%    allLabels - struct with fields
%                   .name: name of synchronization label
%                   .component: components with synchronization label
%                   .location: corresponding locations in each component
%                   .transition: corresponding transition in location
%                   .instant: true/false whether corresponding guard set
%                             yields an instant transition
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Maximilian Perschl, Mark Wetzlinger
% Written:       03-March-2022
% Last update:   19-June-2022 (MW, add warning for lonely syncLabels)
%                23-June-2022 (MW, output argument now struct, not only map)
%                06-April-2023 (MW, rename 'guardempty' -> 'instant')
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init struct for all labels
allLabels = struct([]);

% number of components
numComp = length(pHA.components);

% loop over all components
for i = 1:numComp

    labelsInComponent = {};

    % count in how many components a label occurs
    for j = 1:length(pHA.components(i).location)

        for k = 1:length(pHA.components(i).location(j).transition)
            currentLabel = pHA.components(i).location(j).transition(k).syncLabel;

            if ~isempty(currentLabel)

                % check if corresponding guard set is fullspace -> instant
                % transition
                instant = isa(pHA.components(i).location(j).transition(k).guard,'fullspace');

                % all other entries: check if label already occurred
                found = false;
                if ~isempty(allLabels)
                    found = ismember({allLabels(:).name},currentLabel);
                end
                if any(found)
                    % expand entries corresponding to found name
                    allLabels(found).component = [allLabels(found).component; i];
                    allLabels(found).location = [allLabels(found).location; j];
                    allLabels(found).transition = [allLabels(found).transition; k];
                    allLabels(found).instant = [allLabels(found).instant; instant];
                else
                    % first/new label
                    allLabels = [allLabels; struct('name',currentLabel,'component',i,...
                        'location',j,'transition',k,'instant',instant)];
                end
            end

        end
    end
end

% display warning in case a synchronization label only occurs once (as this
% defeats the purpose of synchronization labels)
onlyOnce = false(length(allLabels),1);
for i=1:length(onlyOnce)
    onlyOnce(i) = length(allLabels(i).component) <= 1;
end
if any(onlyOnce)
    if nnz(onlyOnce) == 1
        warning("Synchronization label '" + allLabels(onlyOnce).name + ...
            "' only occurs once in the entire parallel hybrid automaton.");
    else
        % more than one synchronization label only once
        text = "'" + strjoin({allLabels(onlyOnce).name},"', '") + "'";
        warning("Synchronization labels " + text + ...
            " each only occur once in the entire parallel hybrid automaton.");
    end
end

% ------------------------------ END OF CODE ------------------------------
