function pow = powers(matZ,maxOrder,varargin)
% powers - computes the powers of a matrix zonotope
%
% Syntax:
%    matZ = powers(matZ,maxOrder)
%
% Inputs:
%    matZ - matZonotope object
%    maxOrder - maximum Taylor series order until remainder is computed
%    initialOrder - ?
%    initialPower - ?
%
% Outputs:
%    matZ - matrix zonotope 
%
% Example:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus

% Authors:       Matthias Althoff
% Written:       05-August-2010 
% Last update:   24-September-2010
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% set default values
[initialOrder,initialPower] = setDefaultValues({1,matZ},varargin);

%initialize power
pow{initialOrder}=initialPower;
    
%compute powers
for i=(initialOrder+1):maxOrder
    pow{i} = pow{i-1}*matZ;
end

% ------------------------------ END OF CODE ------------------------------
