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
%    Z = zonotope([0;0],-1+2*rand(2,20));
%    S = dominantDirections(Z);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: 

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

%extract generator matrix
G = Z.G;

%Delete zero-generators
G=nonzeroFilter(G);

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
G=lengthFilter(G,filterLength1);

%apply generator volume filter
Gcells=generatorVolumeFilter(G,filterLength2);

%pick generator with the best volume
Gtemp=volumeFilter(Gcells,Z);
Gpicked=Gtemp{1};

%Select dominant directions S
S(:,1:n)=Gpicked(:,1:n);

% ------------------------------ END OF CODE ------------------------------
