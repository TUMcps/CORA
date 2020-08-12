function [pZ] = project(pZ,dim)
% project - Returns a polyZonotope which is projected onto the specified
%           dimensions
%
% Syntax:  
%    [pZ] = project(pZ,dim)
%
% Inputs:
%    pZ - polyZonotope object
%    dim - projected dimensions
%
% Outputs:
%    pZ - projected polyZonotope
%
% Example: 
%    pZ = polyZonotope([0;0;0],[1 0 0;0 1 0;0 0 1],[0 0;0.1 0;0 0.1],[1 2 3]);
%
%    pZ12 = project(pZ,[1,2]);
%    pZ13 = project(pZ,[1,3]);
%
%    figure
%    plot(pZ12,[1,2],'r','Filled',true,'EdgeColor','none');
%
%    figure
%    plot(pZ13,[1,2],'b','Filled',true,'EdgeColor','none');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/project

% Author:       Niklas Kochdumper
% Written:      25-June-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

pZ.c=pZ.c(dim,:);

if ~isempty(pZ.G)
    pZ.G=pZ.G(dim,:);
end

if ~isempty(pZ.Grest)
    pZ.Grest=pZ.Grest(dim,:);
end

%------------- END OF CODE --------------