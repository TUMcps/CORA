function v = ndimCross(Q)
% ndimCross - Computes the n-dimensional cross product
%
% Syntax:
%    v = ndimCross(Q)
%
% Inputs:
%    Q - matrix of column vectors; must be a n x (n-1) matrix
%
% Outputs:
%    v - n-dimensional cross product 
%
% Example: 
%    Q = rand(4,3);
%    v = ndimCross(Q);
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
for i=1:length(v)
    D=Q;
    D(i,:)=[];
    v(i,1)=(-1)^(i+1)*det(D);
end

% ------------------------------ END OF CODE ------------------------------
