function pow = powers(varargin)
% powers - computes the powers of an interval matrix
%
% Syntax:  
%    [intMat] = powers(varargin) 
%
% Inputs:
%    intMat - interval matrix
%    maxOrder - maximum Taylor series order until remainder is computed
%    initialOrder - first Taylor series order 
%    initialPower - initial power for mixed computations
%
% Outputs:
%    eI - interval matrix exponential
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus

% Author:       Matthias Althoff
% Written:      18-June-2010 
% Last update:  06-July-2010
% Last revision:---

%------------- BEGIN CODE --------------

if nargin==2
    intMat = varargin{1};
    maxOrder = varargin{2};
    initialOrder = 1;
    initialPower = intMat;
elseif nargin==4
    intMat = varargin{1};
    maxOrder = varargin{2};
    initialOrder = varargin{3};
    initialPower = varargin{4};
end

%initialize power
pow{initialOrder}=initialPower;
    
%compute powers
for i=(initialOrder+1):maxOrder
    pow{i} = pow{i-1}*intMat;
end

%------------- END OF CODE --------------