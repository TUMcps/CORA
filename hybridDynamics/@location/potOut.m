function [R] = potOut(obj,R,minInd,maxInd,options)
% potOut - determines the reachable sets after intersection with the
% invariant and obtains the fraction of the reachable set that must have
% transitioned; the resulting reachable sets are all converted to polytopes
%
% Syntax:  
%    [R,endInd] = potOut(obj,R)
%
% Inputs:
%    obj - location object
%    R - cell array of reachable sets
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
%               26-March-2008
%               02-October-2008
%               21-April-2009
%               24-July-2009
%               21-October-2010
%               30-July-2016
%               17-May-2018 (NK, only change sets that intersect guards)
% Last revision:---

%------------- BEGIN CODE --------------

% determine all sets that intersected the guard sets -> sets that are
% partially located outside the invariant
minInd = max(minInd,ones(size(minInd)));
ind = [];
for i = 1:length(minInd)
   temp = minInd(i):maxInd(i);
   ind = [ind,temp];
end
ind = unique(ind);

%parfor (iSet=1:length(R))
for i=1:length(ind)
    
    iSet = ind(i);
    
    if ~iscell(R{iSet})
        
        % overapproximate reachable set by a halfspace representation
        R{iSet} = enclosingPolytope(R{iSet},options);
        
        % intersect with invariant set
        R{iSet} = obj.invariant & R{iSet};      
        
    else
       
        for j = 1:length(R{iSet})
            % overapproximate reachable set by a halfspace representation
            R{iSet}{j} = enclosingPolytope(R{iSet}{j},options);
        
            % intersect with invariant set
            R{iSet}{j} = obj.invariant & R{iSet}{j}; 
        end
    end
end


%------------- END OF CODE --------------