function V = vertices(obj)
% vertices - Computes vertices of an interval object
%
% Syntax:  
%    V = vertices(obj)
%
% Inputs:
%    obj - interval object
%
% Outputs:
%    V - vertices object
%
% Example: 
%    I = interval([1; -1], [2; 1]);
%    V = vertices(I);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/vertices

% Author:       Matthias Althoff
% Written:      24-July-2006 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% compute matrix with all possible generator combinations
I=[1 -1];
for i=1:length(obj.inf)-1
    I=[ones(1,2^i) -ones(1,2^i); I I];
end

Iextended=[ones(1,size(I,2));I];

% convert to zonotope 
zono = zonotope(obj);
Z = zono.Z;

% obtain vertices
V = zeros(size(Z,1),size(I,2));
for i = 1:size(I,2)
    V(:,i) = Z*Iextended(:,i);
end

%------------- END OF CODE --------------