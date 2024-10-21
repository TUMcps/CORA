classdef matrixSet
% matrixSet - abstract superclass for matrix set representations
%
% Syntax:
%    obj = matrixSet()
%
% Inputs:
%    -
%
% Outputs:
%    M - generated matrixSet object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: intervalMatrix, matPolytope, matZonotope

% Authors:       Mark Wetzlinger
% Written:       04-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty abstract class

methods (Access = protected)
    [printOrder] = getPrintSystemInfo(S)
end

end

% ------------------------------ END OF CODE ------------------------------
