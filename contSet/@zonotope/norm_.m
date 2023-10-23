function varargout = norm_(Z,type,mode,varargin)
% norm_ - computes maximum norm value
%
% Syntax:
%    val = norm_(Z,type,mode)
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
%    Z = zonotope([1;0],[1 3 -2; 2 -1 0]);
%    norm(Z)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/norm, minnorm

% Authors:       Victor Gassmann
% Written:       31-July-2020
% Last update:   ---
% Last revision: 27-March-2023 (MW, rename norm_)

% ------------------------------ BEGIN CODE -------------------------------

if strcmp(mode,'exact')
    [varargout{1:2}] = norm_exact(Z,type);
elseif strcmp(mode,'ub')
    varargout{1} = norm_(interval(Z),type);
elseif strcmp(mode,'ub_convex')
    varargout{1} = norm_ub(Z,type);
end

% ------------------------------ END OF CODE ------------------------------
