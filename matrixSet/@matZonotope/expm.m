function eZ = expm(matZono,varargin)
% expm - operator for the exponential matrix of a matrix zonotope
%
% Syntax:  
%    eZ = expm(matZono)
%    eZ = expm(matZono,maxOrder)
%
% Inputs:
%    matZono - matrix zonotope
%    maxOrder - maximum Taylor series order until remainder is computed
%
% Outputs:
%    eZ - matrix zonotope exponential
%
% Example: 
%    C = [0 1;0 -2.5];
%    D = [0 0;0 0.5];
%    intMat = intervalMatrix(C,D);
%    matZono = matZonotope(intMat);
%
%    eZ = expm(matZono)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: intervalMatrix/expm

% Author:       Niklas Kochdumper
% Written:      26-May-2020 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % parse input arguments
    maxOrder = 10;
    
    if nargin > 1
       maxOrder = varargin{1}; 
    end
    
    % compute matrix exponential
    eZ = expmInd(matZono,maxOrder);

%------------- END OF CODE --------------