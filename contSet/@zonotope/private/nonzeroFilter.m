function [G]=nonzeroFilter(G)
% nonZeroFilter - filters out generators of length 0
%
% Syntax:  
%    [G]=nonzeroFilter(G)
%
% Inputs:
%    G - matrix of generators
%
% Outputs:
%    G - reduced matrix of generators
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: 

% Author:        Matthias Althoff, Mark Wetzlinger
% Written:       12-September-2008
% Last update:   02-September-2019
% Last revision: ---

%------------- BEGIN CODE --------------

% delete zero-generators
G = G(:,any(G,1));


%------------- END OF CODE --------------
