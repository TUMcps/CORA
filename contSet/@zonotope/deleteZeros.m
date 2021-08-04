function Z = deleteZeros(Z)
% deleteZeros - removes all generators from a zonotope with zero in all
%    dimensions
%
% Syntax:  
%    Z = deleteZeros(Z)
%
% Inputs:
%    Z - zonotope object
%
% Outputs:
%    Z - zonotope object, every generator has at least one non-zero entry
%
% Example:
%    Z1 = zonotope([0;0],[1 0 -2 0 3 4; 0 0 1 0 -2 1]);
%    Z2 = deleteZeros(Z1);
%    
%    plot(Z1); hold on;
%    plot(Z2,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: 

% Author:        Matthias Althoff
% Written:       15-January-2009
% Last update:   27-Aug-2019
% Last revision: ---

%------------- BEGIN CODE --------------

% extract center and generator matrix
c = center(Z);
G = generators(Z);

% assemble without empty generators
Z.Z = [c,nonzeroFilter(G)];

%------------- END OF CODE --------------
