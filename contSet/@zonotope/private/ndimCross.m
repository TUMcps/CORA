function v = ndimCross(Q)
% ndimCross - Computes the n-dimensional cross product
%
% Syntax:  
%    v = ndimCross(Q)
%
% Inputs:
%    V - matrix of column vectors
%
% Outputs:
%    o - orthogonal vector
%    Q - ???
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

% Author:        Matthias Althoff
% Written:       14-September-2006 
% Last update:   22-March-2007
% Last revision: ---

%------------- BEGIN CODE --------------

v = zeros(length(Q(:,1)),1);
for i=1:length(Q(:,1))
    D=Q;
    D(i,:)=[];
    v(i,1)=(-1)^(i+1)*det(D);
end

%------------- END OF CODE --------------