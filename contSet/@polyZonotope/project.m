function pZ = project(pZ,dims)
% project - projects a polynomial zonotope onto the specified dimensions
%
% Syntax:
%    pZ = project(pZ,dims)
%
% Inputs:
%    pZ - (polyZonotope) polynomial zonotope
%    dims - dimensions for projection
%
% Outputs:
%    pZ - (polyZonotope) projected polynomial zonotope
%
% Example: 
%    pZ = polyZonotope([0;0;0],[1 0 0;0 1 0;0 0 1],[0 0;0.1 0;0 0.1],[1 2 3]);
%
%    pZ_12 = project(pZ,[1,2]);
%    pZ_13 = project(pZ,[1,3]);
%
%    figure;
%    plot(pZ_12,[1,2],'FaceColor','r');
%
%    figure;
%    plot(pZ_13,[1,2],'FaceColor','b');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/project

% Authors:       Niklas Kochdumper
% Written:       25-June-2018
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

pZ.c = pZ.c(dims,:);

if ~isempty(pZ.G)
    pZ.G = pZ.G(dims,:);
end

if ~isempty(pZ.GI)
    pZ.GI = pZ.GI(dims,:);
end

% ------------------------------ END OF CODE ------------------------------
