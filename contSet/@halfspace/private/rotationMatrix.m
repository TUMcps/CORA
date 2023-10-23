function rotMat = rotationMatrix(h, newDir)
% rotationMatrix - computes a rotation matrix to orient the normal vector
% of a hyperplane to newDir
%
% Syntax:
%    rotMat = rotationMatrix(h, newDir)
%
% Inputs:
%    h - halfspace object
%    newDir - vector pointing in the new direction
%
% Outputs:
%    rotMat - rotation matrix
%
% Example: 
%    hs = halfspace([1 -1],2);
%    newDir = [1;1];
%    rotPoint = [0;0];
% 
%    rotate(hs,newDir,rotPoint) % calls rotationMatrix
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus

% Authors:       Matthias Althoff
% Written:       04-September-2013
% Last update:   12-September-2013
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%get dimension
dim_x = length(h.c);
n = h.c/norm(h.c);

if abs(n.'*newDir) ~= 1

    %normalize normal vectors
    newDir = newDir/norm(newDir);
    %create mapping matrix
    B(:,1) = n;
    %find orthonormal basis for n, uVec
    indVec = newDir - (newDir.'*n)*n;
    B(:,2) = indVec/norm(indVec);
    %complete mapping matrix B
    if dim_x>2
        B(:,3:dim_x) = null(B(:,1:2).'); 
    end
    
    %compute angle between uVec and n
    angle = acos(newDir.'*n);
    %rotation matrix
    R = eye(dim_x);
    R(1,1) = cos(angle);
    R(1,2) = -sin(angle);
    R(2,1) = sin(angle);
    R(2,2) = cos(angle);
    %final rotation matrix
    rotMat = B*R*inv(B);
    
else
    if n.'*newDir == 1
        rotMat = eye(dim_x);
    else
        rotMat = -eye(dim_x);
    end
end


% ------------------------------ END OF CODE ------------------------------
