function [totalVol,partialVol] = totalVolumeTest(niP)
% totalVolumeTest - calculates the total volume of the reachable set 
%
% Syntax:
%    [totalVol,partialVol] = totalVolumeTest(niP)
%
% Inputs:
%    niP - non-intersecting polytopes
%
% Outputs:
%    totalVol - joint volume
%    partialVol - volume of each polytope
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
% Last update:   16-August-2007
%                29-October-2007
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%sum volume of all convex polytope parts

%Initialize volume
totalVol=0;
%for each non intersecting polytope
for k=1:length(niP)
    for i=1:length(niP{k})
        partialVol{k}{i}=modVolume(niP{k}{i});
        totalVol=totalVol+partialVol{k}{i};
    end
end

% ------------------------------ END OF CODE ------------------------------
