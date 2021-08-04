function ub = norm_ub(Z,type)
% normbound - computes bound on the maximum norm value
%
% Syntax:  
%    ub = diffnormbound(Z)
%
% Inputs:
%    Z - zonotope object
%
% Outputs:
%   ub      - upper bound on the maximum norm value
%
% Example: 
%   Z = zonotope([randn(2,1),randn(2,10)]);
%   ub = diffnormbound(Z);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: norm

% Author:       Victor Gassmann
% Written:      18-September-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
if ~isYalmipInstalled()
    error('YALMIP must be on the MATLAB search path to use this function');
end
if ~exist('type','var')
    type = 2;
end
if type~=2
    error('Only euclidean norm implemented so far');
end
if ~all(center(Z)==0)
    error('Not implemented yet for non-zero center');
end
%compute upper bound on norm via dual problem of max_{|u|<=1} u'*G'*G*u
G = Z.Z(:,2:end);
[n,m] = size(G);
M = G'*G;
d = sdpvar(m,1);
options = sdpsettings('verbose',0);
%compute optimization problem
optimize([M<=diag(d),diag(d)>=0],d'*ones(m,1),options);
ub = sqrt(value(d'*ones(m,1)));
%------------- END OF CODE --------------