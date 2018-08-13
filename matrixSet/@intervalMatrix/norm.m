function res = norm(obj, varargin)
% norm - computes exactly the maximum norm value of all possible matrices
%
% Syntax:  
%    res = norm(obj, varargin)
%
% Inputs:
%    obj - interval matrix
%    varargin - list of optional inputs
%
% Outputs:
%   res - resulting maximum norm value
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: @zonotope/norm

% Author:       Matthias Althoff
% Written:      02-November-2017
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% convert interval matrix into zonotope matrix
matZ = matZonotope(obj);

% compute norm
res = norm(matZ, varargin{:});

%------------- END OF CODE --------------