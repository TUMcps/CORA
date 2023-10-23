function val = norm(matZ,varargin)
% norm - computes approximately the maximum norm value of all possible
%    matrices
%
% Syntax:
%    val = norm(matZ,varargin)
%
% Inputs:
%    matZ - matZonotope object
%    varargin - list of optional inputs for norm function
%
% Outputs:
%    val - resulting maximum norm value
%
% Example: 
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: @zonotope/norm

% Authors:       Matthias Althoff, Victor Gassmann
% Written:       02-November-2017
% Last update:   23-July-2020
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% convert to interval matrix
M = intervalMatrix(matZ);

% compute norm of interval matrix
val = norm(M,varargin{:});

% ------------------------------ END OF CODE ------------------------------
