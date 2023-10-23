function intMat = intervalMatrix(matP)
% intervalMatrix - computes an enclosing interval matrix of a matrix
%    polytope
%
% Syntax:
%    intMat = intervalMatrix(matP)
%
% Inputs:
%    matP - matPolytope object
%
% Outputs:
%    intMat - intervalMatrix object
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus

% Authors:       Matthias Althoff
% Written:       21-June-2010 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%initialize minimum and maximum values
minMat=Inf(matP.dim);
maxMat=-Inf(matP.dim);

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
C=0.5*(minMat+maxMat);
D=0.5*(maxMat-minMat);
intMat=intervalMatrix(C,D);

% ------------------------------ END OF CODE ------------------------------
