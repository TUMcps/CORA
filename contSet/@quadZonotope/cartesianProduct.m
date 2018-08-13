function qZ = cartesianProduct(qZ1,Z2)
% cartesianProduct - Returns the cartesian product of a quadZonotope and a
% zonotope
%
% Syntax:  
%    qZ = cartesianProduct(qZ1,qZ2)
%
% Inputs:
%    qZ1 - quadZonotope object
%    Z2 - zonotope object
%
% Outputs:
%    qZ - quadZonotope object
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      05-September-2012
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%obtain center and generator of zonotope
Zmat = get(Z2,'Z');
c2 = Zmat(:,1);
G2 = Zmat(:,2:end);
dim2 = length(c2);

%center
c = [qZ1.c; c2];

%first set of generators
G = [qZ1.G; zeros(dim2,length(qZ1.G(1,:)))];

%second set of generators
Gsquare = [];
if ~isempty(qZ1.Gsquare)
    Gsquare = [qZ1.Gsquare; zeros(dim2,length(qZ1.Gsquare(1,:)))];
end

%third set of generators
Gquad = [];
if ~isempty(qZ1.Gquad)
    Gquad = [qZ1.Gquad; zeros(dim2,length(qZ1.Gquad(1,:)))];
end

%fourth set of generators
%determine sizes of generator matrices
[rows_1, cols_1] = size(qZ1.Grest);
[rows_2, cols_2] = size(G2);

%copy generators
Grest(1:rows_1, 1:cols_1) = qZ1.Grest;
Grest((rows_1+1):(rows_1+rows_2), (cols_1+1):(cols_1+cols_2)) = G2;


%generate new zonotope
qZ = quadZonotope(c,G,Gsquare,Gquad,Grest);

%------------- END OF CODE --------------