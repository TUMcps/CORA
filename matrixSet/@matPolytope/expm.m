function eP = expm(matP,varargin)
% expm - operator for the exponential matrix of a matrix polytope
%
% Syntax:
%    eP = expm(matP)
%    eP = expm(matP,maxOrder)
%
% Inputs:
%    matP - matPolytope object
%    maxOrder - maximum Taylor series order until remainder is computed
%
% Outputs:
%    eP - matrix polytope exponential
%
% Example: 
%    C = [0 1;0 -2.5];
%    D = [0 0;0 0.5];
%    intMat = intervalMatrix(C,D);
%    matP = matPolytope(intMat);
%    eP = expm(matP,3)
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
eP = expmInd(matP,maxOrder);

% ------------------------------ END OF CODE ------------------------------
