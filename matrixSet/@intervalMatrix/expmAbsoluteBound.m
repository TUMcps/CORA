function absBound = expmAbsoluteBound(intMat,t)
% expmAbsoluteBound - returns the over-approximation of the absolute bound
% of the symmetric solution of the computation of the mapping mimicing the 
% exponential
%
% Syntax:  
%     absBound = expmAbsoluteBound(intMat,t)
%
% Inputs:
%    intMat - interval matrix
%    t - time
%
% Outputs:
%    absBound - matrix specifying the absolute bound
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus

% Author:       Matthias Althoff
% Written:      02-July-2010 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%extract nominal matrix and symmetric interval matrix
infA = infimum(intMat.int);
supA = supremum(intMat.int);
A = 0.5*(supA+infA);
S = 0.5*(supA-infA);

%set starting order
j=2;

%derive norm values
absA = abs(A);

%generate identity and zero matrices
dim = intMat.dim;
I = eye(dim);
O = zeros(dim);

%build G matrix
G = [absA*t/j O O;...
     S*t/j (absA+S)*t/j O;...
     O I I];

%compute approximated Ginf
%Ginf_test = G^100;

%compute via left and right eigenvectors
%right eigenvector
[Y,eigvalY] = eig(G);
%left eigenvector
[Q,eigvalQ] = eig(G.');
Q = conj(Q);

%select eigenvectors with eigenvalue one
ind_y = find(diag(eigvalY)==1);
ind_q = find(diag(eigvalQ)==1);

%check correct ordering of eigenvectors
%Q.'*Y

%compute Ginf
Ginf = zeros(3*dim);
for i=1:dim
    %get y and q^T
    y = Y(:,ind_y(i));
    qt = Q(:,ind_q(i)).';
    Ginf = Ginf + y*qt/(qt*y);
end

%set initial right hand side
auxMat = [absA*t; S*t; O];

%multiply right hand side with Ginf
auxMat = Ginf*auxMat;

%get absolue bound
absBound = auxMat((2*dim+1):(3*dim),:);

%------------- END OF CODE --------------