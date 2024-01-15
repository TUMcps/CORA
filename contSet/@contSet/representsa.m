function [res,S_conv] = representsa(S,type,varargin)
% representsa - checks if a set can also be represented by a different set,
%    e.g., a special case
%
% Syntax:
%    res = representsa(S,type)
%    res = representsa(S,type,tol)
%    [res,S_conv] = representsa(S,type)
%    [res,S_conv] = representsa(S,type,tol)
%
% Inputs:
%    S - contSet object
%    type - char array
%    tol - (optional) tolerance
%    method - (only conPolyZono) algorithm used for contraction
%             ('forwardBackward', 'linearize', 'polynomial', 'interval', or 'all')
%    iter - (only conPolyZono) number of iteration (integer > 0 or 'fixpoint')
%    splits - (only conPolyZono) number of recursive splits (integer > 0)
%
% Outputs:
%    res - true/false
%    S_conv - converted set
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       19-July-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% default tolerance
[tol,method,splits,iter] = setDefaultValues({1e-12,'linearize',0,1},varargin);

% all admissible comparisons
admissibleTypes = {
    'capsule','conPolyZono','conHyperplane','conZonotope','ellipsoid',...
    'halfspace','interval','levelSet','polytope','polyZonotope',...
    'probZonotope','zonoBundle','zonotope',... % contSet classes
    'origin','point','hyperplane','parallelotope',... % special types
    'emptySet','fullspace'}; % future contSet classes

if isa(S,'conPolyZono')
    % check input arguments
    inputArgsCheck({{S,'att','contSet'};...
                    {type,'str',admissibleTypes};...
                    {tol,'att','numeric',{'scalar','nonnegative'}};
                    {method,'str',{'forwardBackward','linearize',...
                        'polynomial','interval','all'}};
                    {splits,'att','numeric',{'scalar','integer'}};
                    {iter,'att','numeric',{'scalar','integer'}}});

else
    % check input arguments
    inputArgsCheck({{S,'att','contSet'},...
                    {type,'str',admissibleTypes},...
                    {tol,'att','numeric',{'scalar','nonnegative'}}});

end

% call subfunction
if nargout <= 1
    res = representsa_(S,type,tol,method,iter,splits);
elseif nargout == 2
    [res,S_conv] = representsa_(S,type,tol,method,iter,splits);
end

% ------------------------------ END OF CODE ------------------------------
