function res = contains(R,simRes,varargin)
% contains - checks if reachable set contains all simulated
%    points of a set of system trajectory; not for hybrid systems
%
% Syntax:  
%    res = contains(R,simRes)
%    res = contains(R,simRes,type)
%    res = contains(R,simRes,type,tol)
%
% Inputs:
%    R - object of class reachSet
%    simRes - object of class simRes
%    type - (optional) 'exact' or 'approx'
%    tol - (optional) tolerance
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Mark Wetzlinger
% Written:      16-February-2021
% Last update:  28-April-2023 (MW, bug fix for only time-point reachable
%                                  sets, speed up)
% Last revision:---

%------------- BEGIN CODE --------------

% set default values
[type,tol] = setDefaultValues({'exact',1e-10},varargin);

% check input arguments
inputArgsCheck({{R,'att','reachSet'};
                {simRes,'att','simResult'};
                {type,'str',{'exact','approx'}}; ...
                {tol,'att','numeric',{'nonnegative','scalar'}}});

% currently not checked cases
if length(R) > 1
    throw(CORAerror('CORA:notSupported','Multiple branches not supported.'));
end

% check whether time-point or time-interval reachable set given
if ~isempty(R.timeInterval)
    sets = R.timeInterval.set;
    time = R.timeInterval.time;
else
    sets = R.timePoint.set;
    time = R.timePoint.time;
end

% all points have to be checked
nrSimRuns = length(simRes.t);
ptsContained = cell(nrSimRuns,1);
for i=1:nrSimRuns
    ptsContained{i} = false(1,size(simRes.x{i},1));
end

% loop over reachable sets
for k=1:length(sets)
    % for each reachable set, find corresponding simulation points in the
    % entire simResult object
    [pts,ptsContained_] = aux_findPointsInInterval(simRes,time{k});

    % extend list of checked points
    for i=1:nrSimRuns
        ptsContained{i} = ptsContained{i} | ptsContained_{i};
    end

    % exit if a point is found to be outside of the reachable set
    if ~all(contains_(sets{k},pts,type,tol))
        return
    end

end

% check whether all points were checked
if any(cellfun(@(x) ~all(x),ptsContained,'UniformOutput',true))
    throw(CORAerror('CORA:specialError',['Some points of the simulation '...
        'occur outside of the time covered by the reachable set.']));
end

% all good
res = true;

end

% Auxiliary function ------------------------------------------------------

function [pts,ptsContained] = aux_findPointsInInterval(simRes,time)

% tolerance
tol = 1e-12;

% check whether time interval or time point given
ti = isa(time,'interval');

% init points
pts = zeros(size(simRes.x{1},2),0);

% init logical value for contained points
ptsContained = cell(length(simRes.t),1);

% loop over all individual trajectories
for i=1:length(simRes.t)

    % check which points match the correct time
    if ti
        % time points have to be contained in time interval
        ptsContained{i} = contains(time,simRes.t{i}','exact',tol);
    else
        % time points have to be within tolerance of given time point
        ptsContained{i} = withinTol(time,simRes.t{i}',tol);
    end

    % append corresponding state vectors
    pts = [pts simRes.x{i}(ptsContained{i},:)'];
end

end

%------------- END OF CODE --------------