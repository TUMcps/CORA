function eP = expm(matPoly,varargin)
% expm - operator for the exponential matrix of a matrix polytope
%
% Syntax:  
%    eP = expm(matPoly)
%    eP = expm(matPoly,maxOrder)
%
% Inputs:
%    matPoly - matrix polytope
%    maxOrder - maximum Taylor series order until remainder is computed
%
% Outputs:
%    eP - matrix polytope exponential
%
% Example: 
%    C = [0 1;0 -2.5];
%    D = [0 0;0 0.5];
%    intMat = intervalMatrix(C,D);
%    matPoly = matPolytope(intMat);
%
%    eP = expm(matPoly,3)
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
    eP = expmInd(matPoly,maxOrder);

%------------- END OF CODE --------------