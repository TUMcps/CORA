function matI = intervalMatrix(varargin)
% intervalMatrix - computes an enclosing interval matrix of a matrix
% zonotope
%
% Syntax:  
%    matI = intervalMatrix(matZ)
%
% Inputs:
%    matZ - matrix zonotope
%
% Outputs:
%    matI - intervalMatrix
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus

% Author:       Matthias Althoff
% Written:      21-June-2010 
% Last update:  06-October-2010
%               26-August-2011
% Last revision:---

%------------- BEGIN CODE --------------

if nargin==1
    matZ=varargin{1};
    setting=[];
elseif nargin==2
    matZ=varargin{1};
    setting=varargin{2};
end

%delta
delta = abs(matZ.generator{1});
for i=2:matZ.gens
    delta = delta + abs(matZ.generator{i});
end

%instantiate interval matrix
matI = intervalMatrix(matZ.center, delta, setting);


%------------- END OF CODE --------------