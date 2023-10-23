function pow = powers(intMat,varargin)
% powers - computes the powers of an interval matrix
%
% Syntax:
%    pow = powers(intMat,varargin)
%
% Inputs:
%    intMat - intervalMatrix object
%    maxOrder - maximum Taylor series order until remainder is computed
%    initialOrder - first Taylor series order 
%    initialPower - initial power for mixed computations
%
% Outputs:
%    pow - powers of the interval matrix
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus

% Authors:       Matthias Althoff
% Written:       18-June-2010 
% Last update:   06-July-2010
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if nargin==2
    maxOrder = varargin{1};
    initialOrder = 1;
    initialPower = intMat;
elseif nargin==4
    maxOrder = varargin{1};
    initialOrder = varargin{2};
    initialPower = varargin{3};
end

%initialize power
pow{initialOrder}=initialPower;
    
%compute powers
for i=(initialOrder+1):maxOrder
    pow{i} = pow{i-1}*intMat;
end

% ------------------------------ END OF CODE ------------------------------
