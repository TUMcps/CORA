function SpS_out = project(SpS,dims)
% project - projects a spectrahedral shadow onto a set of dimensions
%
% Syntax:
%    SpS_out = project(SpS,dims)
%
% Inputs:
%    SpS - spectraShadow object
%    dims - vector of dimensions
%
% Outputs:
%    SpS_out - projected spectraShadow object
%
% Example: 
%    A0 = eye(4);
%    A1 = [-1 0 0 0;0 1 0 0;0 0 0 0;0 0 0 0];
%    A2 = [0 0 0 0;0 0 0 0;0 0 -1 0;0 0 0 1];
%    SpS = spectraShadow([A0 A1 A2]);
%    SpS_ = project(SpS,1);
%
%    figure; hold on;
%    plot(SpS,[1,2],'Color',[1,0,0,0.5]);
%    plot(SpS_,1,'b','LineWidth',2);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/project

% Authors:       Adrian Kulmburg
% Written:       06-August-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

SpS_out = spectraShadow(SpS.A,SpS.c(dims,:),SpS.G(dims,:));

% Additional properties
SpS_out.emptySet.val = SpS.emptySet.val;
if ~isempty(SpS.bounded.val) && SpS.bounded.val
    SpS_out.bounded.val = true;
end
if ~isempty(SpS.fullDim.val) && SpS.fullDim.val
    SpS_out.fullDim.val = true;
end

if ~isempty(SpS.center.val)
    SpS_out.center.val = SpS.center.val(dims,:);
end

% ------------------------------ END OF CODE ------------------------------
