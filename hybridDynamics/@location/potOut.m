function R = potOut(obj,R,minInd,maxInd,options)
% potOut - determines the reachable sets after intersection with the
% invariant and obtains the fraction of the reachable set that must have
% transitioned; the resulting reachable sets are all converted to polytopes
%
% Syntax:  
%    R = potOut(obj,R,minInd,maxInd,options)
%
% Inputs:
%    obj - location object
%    R - reachSet object storing the reachable set
%    minInd - vector containting the indices of the set which first
%             intersected the guard set for each guard set 
%    maxInd - vector containting the indices of the set which last
%             intersected the guard set for each guard set 
%    options - struct containing algorithm settings
%
% Outputs:
%    R - cell array of reachable sets
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff, Niklas Kochdumper
% Written:      11-May-2007 
% Last update:  18-September-2007
%               21-October-2010
%               30-July-2016
%               17-May-2018 (NK, only change sets that intersect guards)
% Last revision:---

%------------- BEGIN CODE --------------

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
ind = unique(ind);

% loop over all sets that intersect the guard sets
for i = 1:length(ind)
    
    iSet = ind(i);
        
    % overapproximate reachable set by a halfspace representation
    timeInt.set{iSet} = enclosingPolytope(timeInt.set{iSet},options);
    timePoint.set{iSet} = enclosingPolytope(timePoint.set{iSet},options);
        
    % intersect with invariant set
    timeInt.set{iSet} = obj.invariant & timeInt.set{iSet};  
    timePoint.set{iSet} = obj.invariant & timePoint.set{iSet};     
end

% remove last set if it is located outside the invariant
if ~isIntersecting(obj.invariant,timeInt.set{end})
   timeInt.set = timeInt.set(1:end-1); 
   timeInt.time = timeInt.time(1:end-1); 
   timePoint.set = timePoint.set(1:end-1); 
   timePoint.time = timePoint.time(1:end-1); 
end

% construct modified reachSet object
R = reachSet(timePoint,timeInt,R.parent);

%------------- END OF CODE --------------