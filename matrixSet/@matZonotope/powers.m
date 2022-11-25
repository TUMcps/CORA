function pow = powers(varargin)
% powers - computes the powers of a matrix zonotope
%
% Syntax:  
%    matZ = powers(matZ,maxOrder)
%
% Inputs:
%    matZ - matrix zonotope
%    maxOrder - maximum Taylor series order until remainder is computed
%
% Outputs:
%    matZ - matrix zonotope 
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus

% Author:       Matthias Althoff
% Written:      05-August-2010 
% Last update:  24-September-2010
% Last revision:---

%------------- BEGIN CODE --------------

if nargin==2
    matZ = varargin{1};
    maxOrder = varargin{2};
    initialOrder = 1;
    initialPower = matZ;
elseif nargin==4
    matZ = varargin{1};
    maxOrder = varargin{2};
    initialOrder = varargin{3};
    initialPower = varargin{4};
end

%initialize power
pow{initialOrder}=initialPower;
    
%compute powers
for i=(initialOrder+1):maxOrder
    pow{i} = pow{i-1}*matZ;
end


%------------- END OF CODE --------------