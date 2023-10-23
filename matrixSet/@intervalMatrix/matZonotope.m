function matZ = matZonotope(intMat)
% matZonotope - converts an interval matrix to a matrix zonotope
%
% Syntax:
%    matZ = matZonotope(intMat)
%
% Inputs:
%    intMat - intervalMatrix object
%
% Outputs:
%    matZ - matrix zonotope
%
% Example: 
%    c = [2 3 4; 5 6 7];
%    d = [1 0 1; 0 0 1];
%    intMat = intervalMatrix(c,d);
%    matZ = matZonotope(intMat);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff, Mark Wetzlinger
% Written:       21-June-2010 
% Last update:   25-July-2016 (intervalhull replaced by interval)
%                18-June-2023 (MW, fix conversion of interval matrices(!))
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% convert to interval
I = interval(intMat);
% read out dimensions
n = dim(intMat,1);
m = dim(intMat,2);

% read out center
C = center(I);

% compute radius
r = rad(I);

% init generators
G = cell(n*m,1);

% loop over rows and columns
for i=1:n
    for j=1:m
        % linear index
        idx = (i-1)*m+j;
        % init generator with all-zero matrix
        G{idx} = zeros(n,m);
        % overwrite entry using radius
        G{idx}(i,j) = r(i,j);
    end
end
% delete all-zero generators
G(reshape(r'==0,n*m,1)) = [];

% instantiate matrix zonotope
matZ = matZonotope(C,G);

% ------------------------------ END OF CODE ------------------------------
