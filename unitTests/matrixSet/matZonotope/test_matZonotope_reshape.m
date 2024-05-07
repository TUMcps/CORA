function resvec = test_matZonotope_reshape
% test_matZonotope_reshape - unit test function for reshape
% 
% Syntax:
%    res = test_matZonotope_reshape
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Tobias Ladner
% Written:       26-April-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

resvec = [];

% scalar
C = 0;
G = []; G(:,:,1) = 1; G(:,:,2) = -2;
matZscalar = matZonotope(C,G);
dims = {[1,1]};
for i = 1:numel(dims)
    resvec(end+1) = all(dim(reshape(matZscalar,dims{i}(1),dims{i}(2))) == dims{i});
end

% nx1 vector
C = [0; 1; 1];
G = []; G(:,:,1) = [1; -1; -2]; G(:,:,2) = [-2; 0; 1];
matZvector = matZonotope(C,G);
dims = {[1,3],[3,1]};
for i = 1:numel(dims)
    resvec(end+1) = all(dim(reshape(matZvector,dims{i}(1),dims{i}(2))) == dims{i});
end

% matrix
C = [0 2; 1 -1; 1 -2];
G = []; G(:,:,1) = [1 1; -1 0; -2 1]; G(:,:,2) = [-2 0; 0 1; 1 -1];
matZ = matZonotope(C,G);
resvec(end+1) = compareMatrices(C,center(matZ));
dims = {[1,6],[2,3],[6,1]};
for i = 1:numel(dims)
    resvec(end+1) = all(dim(reshape(matZ,dims{i}(1),dims{i}(2))) == dims{i});
end


% combine results
resvec = all(resvec);

% ------------------------------ END OF CODE ------------------------------
