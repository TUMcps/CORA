function P = polytope(matP)
% polytope - Converts a matrix polytope to a polytope 
%
% Syntax:  
%    P = polytope(matP)
%
% Inputs:
%    matP - matrix polytope
%
% Outputs:
%    P - polytope
%
% Example: 
%
% Other m-files required: vertices, polytope
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff
% Written:      21-June-2010
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%convert vertices
for i=1:matP.verts
    V(i,:)=mat2vec(matP.vertex{i});
end
    
%instantiate polytope (using MPT toolbox)
P=polytope(V);




%------------- END OF CODE --------------