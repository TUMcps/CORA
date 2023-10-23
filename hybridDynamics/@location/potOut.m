function R = potOut(loc,R,minInd,maxInd)
% potOut - determines the reachable sets after intersection with the
%    invariant and obtains the fraction of the reachable set that must have
%    transitioned; the resulting reachable sets are all converted to
%    polytopes
%
% Syntax:
%    R = potOut(loc,R,minInd,maxInd)
%
% Inputs:
%    loc - location object
%    R - reachSet object storing the reachable set
%    minInd - vector containing the indices of the set which first
%             intersected the guard set for each guard set 
%    maxInd - vector containing the indices of the set which last
%             intersected the guard set for each guard set
%
% Outputs:
%    R - reachSet object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff, Niklas Kochdumper
% Written:       11-May-2007 
% Last update:   18-September-2007
%                21-October-2010
%                30-July-2016
%                17-May-2018 (NK, only change sets that intersect guards)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% read out time-point and time-interval reachable sets
timeInt = R.timeInterval;
timePoint = R.timePoint;

% determine all sets that intersected the guard sets -> sets that are
% partially located outside the invariant
minInd = max(minInd,ones(size(minInd)));
ind = [];
for i = 1:length(minInd)
    temp = minInd(i):maxInd(i);
    ind = [ind,temp];
end
% remove redundancies
ind = unique(ind);

% loop over all sets that intersect the guard sets
for i=1:length(ind)
    iSet = ind(i);
        
    % overapproximate reachable set by a halfspace representation
    timeInt.set{iSet} = polytope(timeInt.set{iSet});
    timePoint.set{iSet} = polytope(timePoint.set{iSet});
        
    % intersect with invariant set
    timeInt.set{iSet} = and_(loc.invariant,timeInt.set{iSet},'exact');
    timePoint.set{iSet} = and_(loc.invariant,timePoint.set{iSet},'exact');
end

% remove last set if it is located outside the invariant
if ~isIntersecting_(loc.invariant,timeInt.set{end},'exact')
    timeInt.set = timeInt.set(1:end-1); 
    timeInt.time = timeInt.time(1:end-1); 
    timePoint.set = timePoint.set(1:end-1); 
    timePoint.time = timePoint.time(1:end-1);
    % field 'error' currently only supported in linearSys analysis
    if isfield(timeInt,'error')
        timeInt.error = timeInt.error(1:end-1);
    end
    if isfield(timePoint,'error')
        timePoint.error = timePoint.error(1:end-1);
    end
end

% construct modified reachSet object
R = reachSet(timePoint,timeInt,R.parent);

% ------------------------------ END OF CODE ------------------------------
