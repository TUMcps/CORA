function [eI, iPow, E] = expmInd(intMat,varargin)
% expmInd - operator for the exponential matrix of an interval matrix,
%    evaluated independently
%
% Syntax:
%    eI = expmInd(intMat,maxOrder)
%
% Inputs:
%    intMat - interval matrix
%    maxOrder - maximum Taylor series order until remainder is computed
%    initialOrder - ?
%    initialPower - ?
%
% Outputs:
%    eI - interval matrix exponential
%    iPow - cell array storing the powers of the matrix:
%           A,A^2,...,A^(intermediateOrder)
%    E - interval matrix for the remainder
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
%                05-August-2010
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if nargin==2
    
    maxOrder = varargin{1};
    initialOrder = 0;
    initialPower = intMat^0;
    
    %compute powers
    iPow=powers(intMat,maxOrder);
elseif nargin==4
    
    maxOrder = varargin{1};
    initialOrder = varargin{2};
    initialPower = varargin{3};
    
    %compute powers
    iPow=powers(intMat,maxOrder,initialOrder,initialPower);
end

%compute finite Taylor series
%initialize matrix zonotope
eI=initialPower*(1/factorial(initialOrder));
    
%compute Taylor series
for i=(initialOrder+1):maxOrder
    eI = eI + iPow{i}*(1/factorial(i));
end

%compute exponential remainder
E = exponentialRemainder(intMat,maxOrder);

%final result
eI = eI+E;

% ------------------------------ END OF CODE ------------------------------
