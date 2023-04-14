function V = vertices_(I,varargin)
% vertices_ - Computes vertices of an interval object
%
% Syntax:  
%    V = vertices_(I)
%
% Inputs:
%    I - interval object
%
% Outputs:
%    V - vertices object
%
% Example: 
%    I = interval([1; -1], [2; 1]);
%    V = vertices(I);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/vertices_

% Author:       Matthias Althoff
% Written:      24-July-2006 
% Last update:  27-March-2023 (MW, rename vertices_)
% Last revision:05-April-2023 (MW, rewrite to support unbounded intervals)

%------------- BEGIN CODE --------------

% compute all possible combinations of lower/upper bounds
fac = logical(combinator(2,dim(I),'p','r')-1);
nrComb = size(fac,1);

% init all points with lower bound
V = repmat(I.inf,1,nrComb);
% read out supremum
ub = I.sup;

% loop over all factors
for i=1:nrComb
    V(fac(i,:)',i) = ub(fac(i,:));
end

% old version:
% % compute matrix with all possible generator combinations
% one_pm = [1 -1];
% for i=1:length(I.inf)-1
%     one_pm = [ones(1,2^i) -ones(1,2^i); one_pm one_pm];
% end
% 
% Iextended = [ones(1,size(one_pm,2));one_pm];
% 
% % convert to zonotope 
% zono = zonotope(I);
% Z = zono.Z;
% 
% % obtain vertices
% V = zeros(size(Z,1),size(one_pm,2));
% for i = 1:size(one_pm,2)
%     V(:,i) = Z*Iextended(:,i);
% end

%------------- END OF CODE --------------