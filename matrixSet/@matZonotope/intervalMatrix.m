function intMat = intervalMatrix(matZ,varargin)
% intervalMatrix - computes an enclosing interval matrix of a matrix
%    zonotope
%
% Syntax:  
%    intMat = intervalMatrix(matZ)
%    intMat = intervalMatrix(matZ,setting)
%
% Inputs:
%    matZ - matZonotope object
%
% Outputs:
%    intMat - intervalMatrix object
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

% set default values
setting = setDefaultValues({[]},varargin{:});

% center matrix
C = center(matZ);

% delta matrix
D = abs(matZ.generator{1});
for i=2:matZ.gens
    D = D + abs(matZ.generator{i});
end

%instantiate interval matrix
intMat = intervalMatrix(C, D, setting);


%------------- END OF CODE --------------