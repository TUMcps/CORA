function [guards,setIndices] = potInt(loc,R,options)
% potInt - determines which reachable sets potentially intersect with which
%    guard sets
%
% Syntax:  
%    [guards,setIndices] = potInt(loc,R,options)
%
% Inputs:
%    loc - location object
%    R - cell-array of reachable sets
%    options - struct containing the algorithm settings
%
% Outputs:
%    guards - guards that are potentially intersected
%    setIndices - indices of the reachable sets that intersect the guards
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff, Niklas Kochdumper
% Written:      08-May-2007 
% Last update:  26-October-2007
%               20-October-2010
%               27-July-2016
%               23-November-2017
%               03-December-2019 (NK, use approximate intersection test)
% Last revision:---

%------------- BEGIN CODE --------------

% number of reachable sets
nrSets = length(R);
% number of transitions in the location = number of guard sets
nrTrans = length(loc.transition);

% preallocate variables for output arguments (upper bound of entries)
guards = zeros(nrTrans*nrSets,1);
setIndices = zeros(nrTrans*nrSets,1);

% initialize number of intersections
counter = 1;

% loop over all guards
for i = 1:nrTrans
    
    % read out guard set and target location
    guardSet = loc.transition{i}.guard;
    target = loc.transition{i}.target;
    
    % only check if target location is not one of the terminal locations
    if ~all(target == options.finalLoc)
    
        % loop over all reachable sets
        for j = 1:nrSets

            % check if reachable set intersects the guard set
            if isIntersecting_(guardSet,R{j},'approx')
                guards(counter) = i;
                setIndices(counter) = j;
                counter = counter + 1;
            end

        end
    end
end

% remove zeros from resulting lists
guards = guards(1:counter-1);
setIndices = setIndices(1:counter-1);

%------------- END OF CODE --------------