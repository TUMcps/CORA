function cPZ = conPolyZono(C)
% conPolyZono - Converts a capsule to a constrained polynomial zonotope
%
% Syntax:
%    cPZ = conPolyZono(C)
%
% Inputs:
%    C - capsule object
%
% Outputs:
%    cPZ - conPolyZono object
%
% Example: 
%    C = capsule([0;0],[2;2],1);
%    cPZ = conPolyZono(C);
%
%    figure; hold on; xlim([-4,4]); ylim([-4,4]);
%    plot(C,[1,2],'FaceColor','r');
%    plot(cPZ,[1,2],'FaceColor','b','Splits',12);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ellipsoid/conPolyZono

% Authors:       Niklas Kochdumper
% Written:       12-August-2019 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% dimension of capsule
n = dim(C);

% generate constrained polynomial zonotope for the sphere
S = [];
if ~isempty(C.r)
    E = ellipsoid(eye(n)*C.r^2);
    S = conPolyZono(E);
end

% generate constrained polynomial zonotope for the generator
G = [];
if ~isempty(C.g)
    G = conPolyZono(zeros(n,1),C.g,1); 
end

% construct center vector
c = C.c;

if isempty(c)
    c = zeros(n,1); 
end

% combine sphere, generator, and center
if ~representsa(S,'emptySet')
    if ~isempty(G)
        cPZ = G + S + c; 
    else
        cPZ = S + c; 
    end 
else
    if ~isempty(G)
        cPZ = G + s; 
    else
        cPZ = conPolyZono(c,[],[]);
    end
end
    
% ------------------------------ END OF CODE ------------------------------
