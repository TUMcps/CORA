function [totalVol,partialVol]=totalVolumeTest(niP)
% Purpose:  calculates the total volume of the reachable set 
% Pre:      niP (non intersectiong polytopes)
% Post:     total volume
% Tested:   15.09.06,MA
% Modified: 16.08.07,MA
% Modified: 29.10.07,MA

%sum volume of all convex polytope parts---------------
%Initialize volume
totalVol=0;
%for each non intersecting polytope
for k=1:length(niP)
    for i=1:length(niP{k})
        partialVol{k}{i}=modVolume(niP{k}{i});
        totalVol=totalVol+partialVol{k}{i};
    end
end
%------------------------------------------------------