function res = project(cPZ,dims)
% projects - Projects a constrained polynomial zonotope onto the specified
%    dimensions
%
% Syntax:
%    res = project(cPZ,dims)
%
% Inputs:
%    cPZ - (conPolyZono) constrained polynomial zonotope
%    dims - dimensions for projection
%
% Outputs:
%    res - (conPolyZono) projected constrained polynomial zonotope
%
% Example: 
%    c = [0;0;0];
%    G = [2 0 1;0 2 1;1 0 0];
%    E = [1 0 3;0 1 1;0 0 0];
%    A = [1 1 -1.5];
%    b = 0.5;
%    EC = [1 0 0; 0 1 0; 0 0 1];
%    GI = [0.5 0.2; 0 0.1; 0.2 0];
%    cPZ = conPolyZono(c,G,E,A,b,EC,GI);
%     
%    cPZ_ = project(cPZ,[1,2]);
%
%    figure; hold on; box on; grid on;
%    plot(cPZ,[1,2,3],'FaceColor','b','Splits',10);
%    plot(cPZ_,[1,2],'r','Splits',10);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/project, conZonotope/project

% Authors:       Niklas Kochdumper
% Written:       07-November-2018
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check input arguments
% inputArgsCheck({{cPZ,'att','conPolyZono'};
%                 {dims,'att',{'numeric','logical'},{'vector','integer',...
%                     'positive','nonnan'}}});

% project set
res = cPZ;
res.c = cPZ.c(dims);
res.G = cPZ.G(dims,:);

if ~isempty(cPZ.GI)
    res.GI = res.GI(dims,:); 
end
    
% ------------------------------ END OF CODE ------------------------------
