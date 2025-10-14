function S=dominantDirections(Z,varargin)
% dominantDirections - computes the directions that span a parallelotope
%    which tightly encloses a zonotope Z
%
% Syntax:
%    S = dominantDirections(Z,varargin)
%
% Inputs:
%    Z - zonotope object
%    filterLength1 - length of the length filter
%    filterLength2 - length of the generator volume filter
%
% Outputs:
%    S - matrix containing the dominant directions as column vectors
%
% Example:
%    Z = zonotope(zeros(2,1), [ -0.892 -0.097 0.579 0.065 0.743 0.300 -0.848 -0.172 -0.472 0.990 0.562 0.985 -0.152 -0.003 -0.287 0.182 -0.612 0.498 0.893 0.118 ; 0.610 -0.235 -0.271 0.423 -0.343 0.950 0.174 -0.382 0.518 -0.627 -0.608 0.605 0.458 0.618 -0.854 0.820 -0.135 -0.922 0.527 -0.632 ]);
%    S = dominantDirections(Z);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       19-July-2010
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% dimension
n = dim(Z);

% parse input arguments
[filterLength1,filterLength2] = setDefaultValues({n+5,n+3},varargin);

% check input arguments
inputArgsCheck({{Z,'att','zonotope'};
                {filterLength1,'att','numeric',{'nonnan','scalar','positive'}};
                {filterLength2,'att','numeric',{'nonnan','scalar','positive'}}});

%Delete zero-generators
G = nonzeroFilter(Z.G);

% number of generators
nrOfGens = size(G,2);

% if nargin==1
%     filterLength1 = d+5;
%     filterLength2 = d+3;
% elseif nargin==3
%     filterLength1 = varargin{1};
%     filterLength2 = varargin{2};
% end

%correct filter length if necessary
if filterLength1 > nrOfGens
    filterLength1 = nrOfGens;
end

if filterLength2 > nrOfGens
    filterLength2 = nrOfGens;
end

%length filter
G = priv_lengthFilter(G,filterLength1);

%apply generator volume filter
Gcells = priv_generatorVolumeFilter(G,filterLength2);

%pick generator with the best volume
G_picked = priv_volumeFilter(Gcells,Z);

%Select dominant directions S
S(:,1:n) = G_picked{1}(:,1:n);

% ------------------------------ END OF CODE ------------------------------
