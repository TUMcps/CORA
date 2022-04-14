function [tp]=transitionProbabilityTest(niP,field)
% Purpose:  Calculate the transition probability from the actual cell to
%           the reachable cells
% Pre:      niP (non intersecting polytopes), field
% Post:     transition vector
% Tested:   15.09.06,MA
% Modified: 09.10.06,MA


%initialize--------------------------------------------------------
nrOfSegments=get(field,'nrOfSegments');
tp(1:(prod(nrOfSegments)+1),1)=0; %tp: transition probability
%------------------------------------------------------------------

%get total volume of niPs (non intersecting polytopes)-------------
[tv,pv]=totalVolumeTest(niP);
%------------------------------------------------------------------

%get cells that might intersect with the reachable set-------------
for k=1:length(niP)
    for i=1:length(niP{k})
        [tp_total]=cellIntersection2(field,niP{k}{i});
        tp=tp+pv{k}{i}/tv*tp_total;
    end
end