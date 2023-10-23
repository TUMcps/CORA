function tp = transitionProbabilityTest(niP,field)
% transitionProbabilityTest - Calculate the transition probability from the
%    actual cell to the reachable cells
%
% Syntax:
%    tp = transitionProbabilityTest(niP,field)
%
% Inputs:
%    niP - non-intersecting polytopes
%    field - ???
%
% Outputs:
%    tp - transition vector
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Matthias Althoff
% Written:       15-September-2006
% Last update:   09-October-2006
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%initialize
nrOfSegments = get(field,'nrOfSegments');
tp(1:(prod(nrOfSegments)+1),1)=0; %tp: transition probability

%get total volume of niPs (non intersecting polytopes)
[tv,pv] = totalVolumeTest(niP);

%get cells that might intersect with the reachable set
for k=1:length(niP)
    for i=1:length(niP{k})
        [tp_total]=cellIntersection2(field,niP{k}{i});
        tp=tp+pv{k}{i}/tv*tp_total;
    end
end

% ------------------------------ END OF CODE ------------------------------
