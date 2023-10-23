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

% Authors:       Mark Wetzlinger
% Written:       16-February-2021
% Last update:   28-April-2023 (MW, bug fix for only time-point reachable sets, speed up)
%                16-May-2023 (MW, extend to hybrid systems)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% set default values
[type,tol] = setDefaultValues({'exact',1e-10},varargin);

% check input arguments
inputArgsCheck({{R,'att','reachSet'};
                {simRes,'att','simResult'};
                {type,'str',{'exact','approx'}}; ...
                {tol,'att','numeric',{'nonnegative','scalar'}}});

% total number of points: put all points of all trajectories into one big
% list so that contains does not convert the same set into polytopes
% multiple times (costly!)
nrPoints = sum(arrayfun(@(x) sum(cellfun(@length,x.t,'UniformOutput',true)),...
    simRes,'UniformOutput',true));

% logical indexing
ptsContained = false(nrPoints,1);
ptsChecked = false(nrPoints,1);

% loop over locations of reachable set:
% - purely-continuous: location always 0
% - hybrid: different location numbers 1 ... max location
% read out corresponding trajectory parts and sets, check for containment
allLoc = query(R,'allLoc');

for i=1:size(allLoc,2)

    % read out current location
    iLoc = allLoc(:,i);
    
    % read out reachable set of current location
    R_ = find(R,'location',iLoc);

    % loop over all found branches
    for j=1:length(R_)

        % check whether time-point or time-interval reachable set given
        if ~isempty(R_(j).timeInterval)
            sets = R_(j).timeInterval.set;
            time = R_(j).timeInterval.time;
        else
            sets = R_(j).timePoint.set;
            time = R_(j).timePoint.time;
        end
        
        % loop over reachable sets
        for k=1:length(sets)
            % for each reachable set, find corresponding simulation points
            % in the entire simResult object
            [pts,ptsChecked_] = aux_findPointsInInterval(simRes,nrPoints,time{k},iLoc,tol);
        
            % potential improvement: skip points which were already checked
            % and found to be contained

            % check containment
            ptsContained(ptsChecked_) = ptsContained(ptsChecked_) ... 
                | contains_(sets{k},pts,type,tol)';

            % extend list of checked points
            ptsChecked = ptsChecked | ptsChecked_;

            % exit if a point is found to be outside of the reachable set
            % (only if single branch given, otherwise sets from different
            % branches can cover the same time)
            if length(R) == 1 && ~all(ptsContained(ptsChecked_))
                res = false; return
            end
        
        end
    end

    
end

% check whether all points were checked
if ~all(ptsChecked)
    throw(CORAerror('CORA:specialError',['Some points of the simulation '...
        'occur outside of the time covered by the reachable set.']));
end

% check whether all points contained
res = all(ptsContained);

end


% Auxiliary functions -----------------------------------------------------

function [pts,ptsChecked] = aux_findPointsInInterval(simRes,nrPoints,time,loc,tol)

% check whether time interval or time point given
ti = isa(time,'interval');

% init points and logical value for contained points
pts = zeros(size(simRes(1).x{1},2),0);
ptsChecked = false(nrPoints,1);

% index for full list
startIdx = 1;

% loop over all individual trajectories
for r=1:length(simRes)

    % loop over all individual parts
    for part=1:length(simRes(r).t)

        % number of points in this part
        nrPointsPart = length(simRes(r).t{part});

        % skip non-matching locations
        if all(size(simRes(r).loc(part,:)') == size(loc)) ...
                && all(simRes(r).loc(part,:)' == loc)

            % check which points match the correct time
            if ti
                % time points have to be contained in time interval
                tempIdx = contains(time,simRes(r).t{part}','exact',tol);
            else
                % time points have to be within tolerance of given time point
                tempIdx = withinTol(time,simRes(r).t{part}',tol);
            end
            % append checked points
            ptsChecked(startIdx:startIdx+nrPointsPart-1) = tempIdx';
        
            % append corresponding state vectors
            pts = [pts simRes(r).x{part}(tempIdx,:)'];
            
        end

        % shift start index
        startIdx = startIdx + nrPointsPart;

    end
end

end

% ------------------------------ END OF CODE ------------------------------
