function varargout = norm(Z,type,mode)
% norm - computes maximum norm value
%
% Syntax:  
%    val = norm(Z,type,mode)
%
% Inputs:
%    Z    - zonotope object
%    type - (optional) which kind of norm (default: 2)
%    mode - (optional) 'exact', 'ub' (upper bound),'ub_convex' (more
%            precise upper bound computed from a convex program)
%
% Outputs:
%   val - norm value
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: minnorm

% Author:       Victor Gassmann
% Written:      31-July-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
if ~exist('mode','var')
    mode = 'ub';
end
if ~exist('type','var')
    type = 2;
end
if strcmp(mode,'exact')
    [varargout{1:2}] = norm_exact(Z,type);
elseif strcmp(mode,'ub')
    Int = interval(Z);
    varargout{1} = norm(Int,type);
elseif strcmp(mode,'ub_convex')
    varargout{1} = norm_ub(Z,type);
else
    error('Specified type argument not supported');
end
%------------- END OF CODE --------------