function eZ = expm(matZ,varargin)
% expm - operator for the exponential matrix of a matrix zonotope
%
% Syntax:
%    eZ = expm(matZ)
%    eZ = expm(matZ,maxOrder)
%
% Inputs:
%    matZ - matZonotope object
%    maxOrder - maximum Taylor series order until remainder is computed
%
% Outputs:
%    eZ - matrix zonotope exponential
%
% Example: 
%    C = [0;1];
%    D = [0;0.5];
%    intMat = intervalMatrix(C,D);
%    matZ = matZonotope(intMat);
%    eZ = expm(matZ)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: intervalMatrix/expm

% Authors:       Niklas Kochdumper
% Written:       26-May-2020 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input arguments
maxOrder = setDefaultValues({10},varargin);

% compute matrix exponential
eZ = expmInd(matZ,maxOrder);

% ------------------------------ END OF CODE ------------------------------
