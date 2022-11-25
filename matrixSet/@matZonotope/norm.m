function res = norm(obj, varargin)
% norm - computes approximately the maximum norm value of all possible matrices
%
% Syntax:  
%    res = norm(obj, varargin)
%
% Inputs:
%    matZ - matrix zonotope
%    varargin - list of optional inputs
%
% Outputs:
%    res - resulting maximum norm value
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: @zonotope/norm

% Author:       Matthias Althoff, Victor Gassmann
% Written:      02-November-2017
% Last update:  23-July-2020
% Last revision:---

%------------- BEGIN CODE --------------

M = intervalMatrix(obj);
res = norm(M,varargin{:});

%------------- END OF CODE --------------