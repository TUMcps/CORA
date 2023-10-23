function R = order(R)
% order - orders the elements of a reachSet object chronologically;
%    currently only supported for one branch
%
% Syntax:
%    R = order(R)
%
% Inputs:
%    R - reachSet object
%
% Outputs:
%    R - reachSet object
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       07-November-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% only supported for one branch
if length(R) > 1
    throw(CORAerror('CORA:notSupported','Only supports one branch.'));
end

% time-point solutions
time = cell2mat(R.timePoint.time);
[time,idx] = sort(time,'ascend');
R.timePoint.set = R.timePoint.set(idx);
R.timePoint.time = num2cell(time);
% error only for linearSys
if isfield(R.timePoint,'error')
    R.timePoint.error = R.timePoint.error(idx);
end

% time-interval solutions
if ~isempty(R.timeInterval)
    time = R.timeInterval.time;
    nrSets = length(time);
    % read out interval
    timeInf = zeros(nrSets,1);
    timeSup = zeros(nrSets,1);
    for i=1:nrSets
        timeInf(i,1) = time{i}.infimum;
        timeSup(i,1) = time{i}.supremum;
    end
    % order depending on start time of interval
    [timeInf,idx] = sort(timeInf,'ascend');
    timeSup = timeSup(idx);
    
    % new list of sets
    R.timeInterval.set = R.timeInterval.set(idx);
    % new list of time intervals
    for i=1:nrSets
        R.timeInterval.time{i} = interval(timeInf(i),timeSup(i));
    end
    % error only for linearSys
    if isfield(R.timeInterval,'error')
        R.timeInterval.error = R.timeInterval.error(idx);
    end
end

% ------------------------------ END OF CODE ------------------------------
