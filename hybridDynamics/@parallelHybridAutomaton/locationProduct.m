function loc = locationProduct(pHA,locID,varargin)
% locationProduct - constructs an overall location object from the active 
%    locations of the subcomponents with a local automaton product
%
% Syntax:
%    loc = locationProduct(pHA,locID)
%    loc = locationProduct(pHA,locID,allLabels)
%
% Inputs:
%    pHA - parallelHybridAutomaton object
%    locID - IDs of the currently active locations
%    allLabels - information about synchronization labels
%
% Outputs:
%    loc - constructed location object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Johann Schoepfer, Niklas Kochdumper
% Written:       08-June-2018  
% Last update:   09-July-2018
%                10-October-2024 (MW, use only pHA and locID for merge)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% information about labels may be given or not
narginchk(2,3);

% merge transition sets
[mergedTransSets,mergedLabels] = mergeTransitionSets(pHA,locID,varargin{:});

% merge invariants
mergedInvSet = mergeInvariants(pHA,locID,mergedLabels);

% merge flows
mergedFlow = mergeFlows(pHA,locID);

% construct resulting location object
loc = location(mergedInvSet,mergedTransSets,mergedFlow);

% ------------------------------ END OF CODE ------------------------------
