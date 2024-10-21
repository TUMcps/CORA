function SpS = spectraShadow(zB)
% spectraShadow - converts a zonotope bundle to a spectrahedral shadow
%
% Syntax:
%    SpS = spectraShadow(zB)
%
% Inputs:
%    zB - zonoBundle object
%
% Outputs:
%    SpS - spectraShadow object
%
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Adrian Kulmburg
% Written:       19-August-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if zB.parallelSets == 0
    SpS = spectraShadow(); return
end

% transform all zonotopes to spectrahedral shadows
SpStmp = cell(zB.parallelSets,1);
for i=1:zB.parallelSets
    SpStmp{i} = spectraShadow(zB.Z{i});
end

% intersect interval hulls
SpS = SpStmp{1};
for i=2:zB.parallelSets
    SpS = and_(SpS,SpStmp{i},'exact');
end

% Additional properties
SpS.bounded.val = true;
SpS.emptySet.val = representsa_(zB,'emptySet',1e-10);
SpS.fullDim.val = isFullDim(zB);
SpS.center.val = center(zB);

% ------------------------------ END OF CODE ------------------------------
