function varargout = norm(Z,varargin)
% norm - computes maximum norm value
%
% Syntax:  
%    val = norm(Z,type,mode)
%
% Inputs:
%    Z - zonotope object
%    type - (optional) which kind of norm (default: 2)
%    mode - (optional) 'exact', 'ub' (upper bound),'ub_convex' (more
%            precise upper bound computed from a convex program)
%
% Outputs:
%    val - norm value
%    x - vertex attaining maximum norm
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

% pre-processing
[res,vars] = pre_norm('zonotope',Z,varargin{:});

% check premature exit
if res
    % if result has been found, it is stored in the first entry of var
    varargout{1} = vars{1}; return
else
    % assign values
    Z = vars{1}; type = vars{2}; mode = vars{3};
end


if strcmp(mode,'exact')
    [varargout{1:2}] = norm_exact(Z,type);
elseif strcmp(mode,'ub')
    varargout{1} = norm(interval(Z),type);
elseif strcmp(mode,'ub_convex')
    varargout{1} = norm_ub(Z,type);
end

%------------- END OF CODE --------------