function Z = generateRandom(varargin)
% generateRandom - Generates a random zonotope
%
% Syntax:  
%    Z = generateRandom(varargin)
%
% Inputs:
%    dim      - (optional) dimension
%    cen      - (optional) center
%    nrOfGens - (optional) number of generators
%    type     - (optional) distribution from which generators are created
%                'uniform', 'exp', 'gamma'
%
% Outputs:
%    Z - random zonotope
%
% Example: 
%    Z1 = zonotope.generateRandom();
%    Z2 = zonotope.generateRandom(3);
%    Z3 = zonotope.generateRandom(2,ones(2,1));
%    Z4 = zonotope.generateRandom(4,zeros(4,1),10);
%    Z5 = zonotope.generateRandom(3,zeros(3,1),20,'gamma');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Matthias Althoff, Mark Wetzlinger
% Written:      17-Sep-2019
% Last update:  24-Sep-2019
%               01-May-2020 (integration of randomZonotope.m)
% Last revision:---

%------------- BEGIN CODE --------------

% define bounds for random values
dim_low = 2;
dim_up = 10;
orderGen_low = 1;
orderGen_up = 10;
center_low = -5;
center_up = 5;
maxLength = 5; % only type = 'uniform'

if nargin < 4
    type = 'uniform';
else
    type = varargin{4};
end

% generate random values
if nargin >= 1 && ~isempty(varargin{1})
    n = varargin{1};
else
    n = dim_low + floor(rand(1) * (dim_up - dim_low + 1));
end

c = center_low + rand(n,1) * (center_up - center_low);

numGen_low = floor(orderGen_low*n);
numGen_up = floor(orderGen_up*n);

nrOfGens = numGen_low + floor(rand(1) * (numGen_up - numGen_low));

if nargin >= 2 && ~isempty(varargin{2})
   c = varargin{2};
end
if nargin >= 3 && ~isempty(varargin{3})
   nrOfGens = varargin{3}; 
end

%generate random vector for the generator lengths
l = zeros(nrOfGens,1);
%uniform distribution
if strcmp(type,'uniform')
    l = maxLength*rand(nrOfGens,1);

%exponential distribution    
elseif strcmp(type,'exp')
    l = exprnd(1,nrOfGens,1);

%gaussian distribution    
elseif strcmp(type,'gamma')
    l = gamrnd(2,1,nrOfGens,1);
end

% init generator matrix
G = zeros(n,nrOfGens);
%create generators
for i=1:nrOfGens
    %generate random point on sphere
    gTmp = randomPointOnSphere(n);
    %stretch by length
    G(:,i) = l(i)*gTmp;
end

% instantiate zonotope
Z = zonotope(c,G);

end

%------------- END OF CODE --------------