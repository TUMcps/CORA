function Z = unify(Z1, Z2)
% or - Computes union of zonotopes; works best when directions of
% generators are similar
%
% Syntax:  
%    Z = or(Z1, Z2)
%
% Inputs:
%    Z1 - zonotope object
%    Zcell - cell array of zonotope objects
%
% Outputs:
%   Z - zonotope object
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff
% Written:      20-September-2013
% Last update:  25-July-2016 (intervalhull replaced by interval)
% Last revision:---

%------------- BEGIN CODE --------------


%get Z-matrix from zonotope Z
Zmatrix1=get(Z1,'Z');
Zmatrix2=get(Z2,'Z');

%extract generator matrix
G1=Zmatrix1(:,2:end);
G2=Zmatrix2(:,2:end);

%obtain matrix of points from generator matrix
V = [G1,-G1, G2, -G2];

%compute the arithmetic mean of the vertices
mean=sum(V,2)/length(V(1,:));

%obtain sampling matrix
translation=mean*ones(1,length(V(1,:)));
sampleMatrix=V-translation;

%compute the covariance matrix
C=cov(sampleMatrix');

%singular value decomposition
[U,S,V] = svd(C);

%auxiliary computations
orientedMatrix=U'*sampleMatrix;

%Zred1
Ztrans = U.'*Z1;
Zinterval1 = interval(Ztrans);

%Zred2
Ztrans = U.'*Z2;
Zinterval2 = interval(Ztrans);

%Zinterval
Zinterval = Zinterval1 | Zinterval2;

%enclosure
Z = U*zonotope(Zinterval);


%------------- END OF CODE --------------