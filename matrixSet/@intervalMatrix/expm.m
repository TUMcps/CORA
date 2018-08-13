function [eI, iPow, E] = expm(varargin)
% expm - operator for the exponential matrix of an 
% interval matrix, evaluated dependently
%
% Syntax:  
%    eI = expmInd(intMat,maxOrder)
%
% Inputs:
%    intMat - interval matrix
%    maxOrder - maximum Taylor series order until remainder is computed
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
%               05-August-2010
% Last revision:---

%------------- BEGIN CODE --------------

if nargin==3
    intMat = varargin{1};
    r = varargin{2};
    maxOrder = varargin{3};
    
    %compute exact terms
    [sq,H] = dependentTerms(intMat*(1/r),r);
    
    initialOrder = 2;
    initialPower = sq;
    
    %init eI
    eI = H;
    
elseif nargin==5
    intMat = varargin{1};
    r = varargin{2};
    maxOrder = varargin{3};
    initialOrder = varargin{4};
    initialPower = varargin{5};
    
    %init eI
    eI=initialPower*(1/factorial(initialOrder));
end

%compute powers
iPow=powers(intMat,maxOrder,initialOrder,initialPower);
    
%compute Taylor series
for i=(initialOrder+1):maxOrder
    eI = eI + iPow{i}*(1/factorial(i));
end

%compute exponential remainder
E = exponentialRemainder(intMat,maxOrder);

%final result
eI = eI+E;

%------------- END OF CODE --------------