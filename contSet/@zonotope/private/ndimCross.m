function v = ndimCross(Q)
% ndimCross - Computes the n-dimensional cross product
%
% Syntax:
%    v = ndimCross(Q)
%
% Inputs:
%    Q - matrix of column vectors
%
% Outputs:
%    v - n-dimensional cross product 
%
% Example: 
%    Z=zonotope(rand(2,5));
%    P=polytope(Z);
%    plot(P);
%    hold on
%    plot(Z);
%
% Other m-files required: vertices, polytope
% Subfunctions: none
% MAT-files required: none
%
% See also: interval,  vertices

% Authors:       Matthias Althoff
% Written:       14-September-2006 
% Last update:   22-March-2007
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

dim = length(Q(:,1));
v = zeros(dim,1);
for i=1:dim
    D=Q;
    D(i,:)=[];
    v(i,1)=(-1)^(i+1)*det(D);
end

% ------------------------------ END OF CODE ------------------------------
