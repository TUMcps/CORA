function matI = intervalMatrix(matP)
% intervalMatrix - computes an enclosing interval matrix of a matrix
% polytope
%
% Syntax:  
%    matI = intervalMatrix(matP)
%
% Inputs:
%    matP - matrix polytope
%
% Outputs:
%    matI - intervalMatrix
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus

% Author:       Matthias Althoff
% Written:      21-June-2010 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%initialize minimum and maximum values
minMat=ones(matP.dim)*inf;
maxMat=ones(matP.dim)*(-inf);

%find minimum and maximum values
for i=1:matP.verts
    %find smaller values
    ind=find(matP.vertex{i}<minMat);
    %update minimum values
    minMat(ind)=matP.vertex{i}(ind);
    
    %find greater values
    ind=find(matP.vertex{i}>maxMat);
    %update minimum values
    maxMat(ind)=matP.vertex{i}(ind);   
end

%instantiate interval matrix
center=0.5*(minMat+maxMat);
delta=0.5*(maxMat-minMat);
matI=intervalMatrix(center,delta);

%------------- END OF CODE --------------