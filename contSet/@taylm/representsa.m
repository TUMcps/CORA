function [res,t_conv] = representsa(t,type,varargin)
% representsa - checks if a set can also be represented by a different set,
%    e.g., a special case
%
% Syntax:
%    res = representsa(t,type)
%    res = representsa(t,type,tol)
%    [res,t_conv] = representsa(t,type)
%    [res,t_conv] = representsa(t,type,tol)
%
% Inputs:
%    t - taylm object
%    type - char array
%    tol - (optional) tolerance
%    method - (only conPolyZono) algorithm used for contraction
%             ('forwardBackward', 'linearize', 'polynomial', 'interval', or 'all')
%    iter - (only conPolyZono) number of iteration (integer > 0 or 'fixpoint')
%    splits - (only conPolyZono) number of recursive splits (integer > 0)
%
% Outputs:
%    res - true/false
%    t_conv - converted set
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Adrian Kulmburg
% Written:       19-January-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% default tolerance
[tol,method,splits,iter] = setDefaultValues({1e-12,'linearize',0,1},varargin);

% all admissible comparisons
admissibleTypes = {
    'capsule','conPolyZono','conHyperplane','conZonotope','ellipsoid',...
    'halfspace','interval','levelSet','polygon','polytope','polyZonotope',...
    'probZonotope','zonoBundle','zonotope',... % contSet classes
    'origin','point','hyperplane','parallelotope','convexSet',... % special types
    'emptySet','fullspace'};

% check input arguments
inputArgsCheck({{t,'att','taylm'},...
    {type,'str',admissibleTypes},...
    {tol,'att','numeric',{'scalar','nonnegative'}}});


% call subfunction, by first converting to a polynomial zonotope
pZ = polyZonotope(t);
if nargout <= 1
    res = representsa_(pZ,type,tol,method,iter,splits);
elseif nargout == 2
    [res,pZ_conv] = representsa_(pZ,type,tol,method,iter,splits);
    % convert back to taylm
    t_conv = taylm(pZ_conv);
end

% ------------------------------ END OF CODE ------------------------------
