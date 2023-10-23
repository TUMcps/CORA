function R = initReachSet(timePoint, timeInt)
% initReachSet - create and object of class reachSet that stores the
%    reachable set
%
% Syntax:
%    R = reachSet.initReachSet(timePoint,timeInt)
%
% Inputs:
%    timePoint - struct containing the time-point reachable set
%    timeInt - struct containing the time-interval reachable set
%
% Outputs:
%    R - reachSet object
%
% Example:
%    -
%
% References:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: previously helper function 'createReachSetObject'

% Authors:       Niklas Kochdumper
% Written:       ---
% Last update:   19-November-2022 (MW, move from 'createReachSetObject')
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% no time-interval solution (e.g., observer)
if nargin == 1
    timeInt = [];
end

% remove empty cells (occurs due to instantiation of full list and
% premature exit, e.g., because of specification violation)
checkEmpty = @(x) ~isa(x,'contSet') && isempty(x);
if ~isempty(timeInt) && checkEmpty(timeInt.set{end})
    ind = ~cellfun(@(x) checkEmpty(x), timeInt.set);
    timeInt.set = timeInt.set(ind);
    timeInt.time = timeInt.time(ind);
    if isfield(timeInt, 'algebraic')
        timeInt.algebraic = timeInt.algebraic(ind);
    end
    if isfield(timeInt, 'error')
        timeInt.error = timeInt.error(ind);
    end
end
if ~isempty(timePoint) && checkEmpty(timePoint.set{end})
    ind = ~cellfun(@(x) checkEmpty(x), timePoint.set);
    timePoint.set = timePoint.set(ind);
    timePoint.time = timePoint.time(ind);
    if isfield(timePoint, 'error')
        timePoint.error = timePoint.error(ind);
    end
end

% check if sets are split: then, timeInt.set{1} would be a cell-array of
% all parallel sets (always the case for nonlinear system classes)
splitStructure = (~isempty(timePoint) && iscell(timePoint.set{1})) ...
    || (~isempty(timeInt) && iscell(timeInt.set{1}));

if ~splitStructure
    if isa(timePoint.set{1}, 'zonoBundle') || ~representsa_(timePoint.set{1}, 'emptySet', eps)
        R = reachSet(timePoint, timeInt);
    else
        R = reachSet([], timeInt);
    end
    return;
end

% nonlinear systems have the structure of split sets, that is,
%   timePoint.set{i}{j}, (or timeInt.set{i}{j})
% where i is the time step and j is the j-th parallel set, but j can be 1,
% and the sets are not actually split; thus, we check whether there are
% multiple sets at the end to check if splitting occurred at any point
noSplitOccurred = (~isempty(timePoint) && length(timePoint.set{end}) == 1) ...
    || (~isempty(timeInt) && length(timeInt.set{end}) == 1);

if noSplitOccurred
    % no splitting occurred -> copy sets and exit
    if ~isempty(timePoint.set)
        timePoint.set = cellfun(@(x) x{1}.set, timePoint.set, 'UniformOutput', false);
    end
    if ~isempty(timeInt.set)
        timeInt.set = cellfun(@(x) x{1}, timeInt.set, 'UniformOutput', false);
    end
    if isfield(timeInt, 'algebraic')
        timeInt.algebraic = cellfun(@(x) x{1}, timeInt.algebraic, 'UniformOutput', false);
    end
    R = reachSet(timePoint, timeInt);
    return;
end

% split sets -> bring to correct format
timeInt_ = {};
timePoint_ = {};

ind = 1:length(timePoint.set{2});
parent = zeros(length(timePoint.set{2}), 1);

% loop over all time steps
for i = 1:length(timeInt.set)
    % i is shifted by +1 for timePoint

    % modify index vector
    ind_ = zeros(length(timePoint.set{i+1}), 1);
    maxInd = max(ind);

    for j = 1:length(timePoint.set{i+1})
        if isfield(timePoint.set{i+1}{j}, 'parent') && i > 1
            ind_(j) = maxInd + 1;
            parent = [parent; ind(timePoint.set{i+1}{j}.parent)];
            maxInd = maxInd + 1;
        else
            ind_(j) = ind(timePoint.set{i+1}{j}.prev);
        end
    end
    ind = ind_;

    % handling of initial set
    if i == 1
        for j = 1:max(ind)
            timePoint_{j}.set = {timePoint.set{1}{1}.set};
            timePoint_{j}.time = timePoint.time(1);
        end
    end

    % copy entries
    for j = 1:length(timeInt.set{i})
        if ind(j) <= length(timeInt_)
            % append to end of branch number 'ind'
            timeInt_{ind(j)}.set = [timeInt_{ind(j)}.set; timeInt.set{i}(j)];
            timeInt_{ind(j)}.time = [timeInt_{ind(j)}.time; timeInt.time(i)];
            if isfield(timeInt, 'algebraic')
                timeInt_{ind(j)}.algebraic = [timeInt_{ind(j)}.algebraic; timeInt.algebraic{i}(j)];
            end
        else
            % extend cell-array by another branch: number 'ind(j)'
            timeInt_{ind(j)}.set = timeInt.set{i}(j);
            timeInt_{ind(j)}.time = timeInt.time(i);
            if isfield(timeInt, 'algebraic')
                timeInt_{ind(j)}.algebraic = timeInt.algebraic{i}(j);
            end
        end
        if ind(j) <= length(timePoint_)
            % append to end of branch number 'ind'
            timePoint_{ind(j)}.set = [timePoint_{ind(j)}.set; {timePoint.set{i+1}{j}.set}];
            timePoint_{ind(j)}.time = [timePoint_{ind(j)}.time; timePoint.time(i+1)];
        else
            % extend cell-array by another branch: number 'ind(j)'
            timePoint_{ind(j)}.set = {timePoint.set{i+1}{j}.set};
            timePoint_{ind(j)}.time = timePoint.time(i+1);
        end
    end
end

% generate reachSet object
for i = 1:length(timeInt_)

    R_ = reachSet(timePoint_{i}, timeInt_{i});

    if i == 1
        R = R_;
    else
        if parent(i) > 0
            R = add(R, R_, parent(i));
        else
            R = add(R, R_);
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
