function Gred = priv_kmeansFilter(G,rem)
% priv_kmeansFilter - filters out generators by the k-means clustering
%    algorithm
%
% Syntax:
%    Gred = priv_kmeansFilter(G,rem)
%
% Inputs:
%    G - matrix of generators
%    rem - number of remaining generators
%
% Outputs:
%    Gred - reduced matrix of generators
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       12-September-2008
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% filter generators with k-means
[~,C] = kmeans(G',rem);
Gred = C';

% ------------------------------ END OF CODE ------------------------------
