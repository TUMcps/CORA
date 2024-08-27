function v = ndimCross(Q)
% ndimCross - computes the n-dimensional cross product acc. to [1, Sec. 5.1]
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
%    Q = [1 2; 3 4; 5 6];
%    v = ndimCross(Q);
%
% Reference:
%    [1] M. Althoff, O. Stursberg, M. Buss. "Computing reachable sets of
%        hybrid systems using a combination of zonotopes and polytopes",
%        Nonlinear Analysis: Hybrid Systems 4 (2010), p. 233-249.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/polytope

% Authors:       Matthias Althoff, Mark Wetzlinger
% Written:       14-September-2006 
% Last update:   22-March-2007
%                01-August-2024 (MW, simplify computation)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

v = zeros(size(Q,1),1);
for i=1:length(v)
    v(i,1) = (-1)^(i+1) * det(Q([1:i-1,i+1:end],:));
end

% ------------------------------ END OF CODE ------------------------------
