function pHA = priv_derivatives(pHA,locID)
% priv_derivatives - compute derivatives of reset functions for all 
%    transitions of a given location product object (must be first
%    computed!)
%    note: we currently compute the derivatives in mergeTransitionSets
%
% Syntax:
%    pHA = priv_derivatives(pHA,locID)
%
% Inputs:
%    pHA - parallelHybridAutomaton object
%    locID - location ID vector
%
% Outputs:
%    pHA - updated parallelHybridAutomaton object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: location/derivatives, transition/derivatives,
%    nonlinearReset/derivatives

% Authors:       Mark Wetzlinger
% Written:       15-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% private function -- no input arguments check

% find location with correct identifier in pHA object
idxInLocProd = find(all([pHA.locProd.locID] == locID,1),1);
if isempty(idxInLocProd)
    throw(CORAerror('CORA:specialError',...
        'No such location product has been computed.'));
end
% (note: idxInLocProd should have exactly length 1)
loc = pHA.locProd(idxInLocProd(1));

% all transitions
transIdx = 1:numel(loc.transition);
% set correct path
fpath = [CORAROOT filesep 'models' filesep 'auxiliary' filesep pHA.name ...
    filesep 'location_' replace(num2str(reshape(locID,1,[])),' ','')];
% call location/derivatives
pHA.locProd.location(idxInLocProd) = derivatives(loc,transIdx,fpath);

% ------------------------------ END OF CODE ------------------------------
