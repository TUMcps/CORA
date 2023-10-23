function Z = exactPlus(Z,Z2,varargin)
% exactPlus - adds two zonotopes by adding all generators of common
%    generator factors. Caution: The correspondance of generator factors
%    has to be ensured before calling the function; this function is not a
%    replacement for the Minkowski sum
%
% Syntax:
%    Z = exactPlus(Z,Z2)
%
% Inputs:
%    Z - zonotope object
%    Z2 - zonotope object
%    nrOfGens - (optional) limit on the number of generators that can be
%               added exactly
%
% Outputs:
%    Z - zonotope object
%
% Example: 
%    Z1 = zonotope([0;0],[1 2 -1; 1 -1 3]);
%    Z2 = zonotope([0;0],[2 4 -1; 2 -2 3]);
%    Zexact = exactPlus(Z1,Z2);
%    Z = Z1 + Z2;
% 
%    figure; hold on;
%    plot(Z1);
%    plot(Z2);
%    plot(Zexact,[1,2],'g');
%    plot(Z,[1,2],'r--');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/plus

% Authors:       Matthias Althoff
% Written:       30-August-2013
% Last update:   06-September-2013
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% number of generators
nrOfgens1 = size(Z.G,2);
nrOfgems2 = size(Z2.G,2);

if nargin == 2
    nrOfGens = min([nrOfgens1, nrOfgems2]);
elseif nargin == 3
    nrOfGens = min([nrOfgens1, nrOfgems2, varargin{1}]);
end

% resulting zonotope
Z.c = Z.c + Z2.c;
Z.G = [Z.G(:,1:nrOfGens) + Z2.G(:,1:nrOfGens),...
    Z.G(:,nrOfGens+1:end), ...
    Z2.G(:,nrOfGens+1:end)];

% ------------------------------ END OF CODE ------------------------------
