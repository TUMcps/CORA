function S=dominantDirections(varargin)
% dominantDirections - computes the directions that span a parallelotope
% which tightly encloses a zonotope Z
%
% Syntax:  
%    S=dominantDirections(varargin)
%
% Inputs:
%    Z - zonotope object
%    filterLength1 - length of the length filter
%    filterLength2 - length of the generator volume filter
%
% Outputs:
%    S - matrix containing the dominant directions as column vectors
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: 

% Author:       Matthias Althoff
% Written:      19-July-2010
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%get Z-matrix from zonotope Z
Z = varargin{1};
Zmatrix=get(Z,'Z');
dim=length(Zmatrix(:,1));

%extract generator matrix
G=Zmatrix(:,2:end);

%Delete zero-generators
G=nonzeroFilter(G);

if nargin==1
    filterLength1 = dim+5;
    filterLength2 = dim+3;
elseif nargin==3
    filterLength1 = varargin{2};
    filterLength2 = varargin{3};
end

%correct filter length if necessary
if filterLength1>length(G(1,:))
    filterLength1=length(G(1,:));
end

if filterLength2>length(G(1,:))
    filterLength2=length(G(1,:));
end

%length filter
G=lengthFilter(G,filterLength1);

%apply generator volume filter
Gcells=generatorVolumeFilter(G,filterLength2);

%pick generator with the best volume
Gtemp=volumeFilter(Gcells,Z);
Gpicked=Gtemp{1};

%Select dominant directions S
S(:,1:dim)=Gpicked(:,1:dim);



%------------- END OF CODE --------------
