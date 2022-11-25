function [eI, iPow, E] = expm(varargin)
% expm - operator for the exponential matrix of an 
% interval matrix, evaluated dependently
%
% Syntax:  
%    eI = expm(intMat)
%    eI = expm(intMat,maxOrder)
%    [eI, iPow, E] = expm(intMat,r,maxOrder)
%
% Inputs:
%    intMat - interval matrix
%    maxOrder - maximum Taylor series order until remainder is computed
%    r - time step size
%
% Outputs:
%    eI - interval matrix exponential
%    iPow - cell array storing the powers of the matrix:
%           A,A^2,...,A^(intermediateOrder)
%    E - interval matrix for the remainder
%
% Example: 
%    C = [0 1;0 -2.5];
%    D = [0 0;0 0.5];
%    intMat = intervalMatrix(C,D);
%
%    eI = expm(intMat)
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

if nargin == 1
    
    intMat = varargin{1};
    maxOrder = 10;
    
    %compute exact terms
    [sq,H] = dependentTerms(intMat,1);
    
    initialOrder = 2;
    initialPower = sq;
    
    %init eI
    eI = H;

elseif nargin == 2
    
    intMat = varargin{1};
    maxOrder = varargin{2};
    
    %compute exact terms
    [sq,H] = dependentTerms(intMat,1);
    
    initialOrder = 2;
    initialPower = sq;
    
    %init eI
    eI = H;
    

elseif nargin==3
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