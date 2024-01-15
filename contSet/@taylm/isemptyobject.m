function res = isemptyobject(tay)
% isemptyobject - checks whether a Taylor model contains any information at
%    all; consequently, the set is interpreted as the empty set 
%
% Syntax:
%    res = isemptyobject(tay)
%
% Inputs:
%    tay - taylm object
%
% Outputs:
%    res - true/false
%
% Example: 
%    tay = taylm(interval(1,2),4,'x');
%    isemptyobject(tay); % false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       24-July-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = isnumeric(tay.coefficients) && isscalar(tay.coefficients) && tay.coefficients == 0 ...
        && isnumeric(tay.monomials) && isempty(tay.monomials) ...
        && isa(tay.remainder,'interval') ...
        && isscalar(tay.remainder.inf) && isscalar(tay.remainder.sup) ...
        && tay.remainder.inf == 1 && tay.remainder.sup == 1 ...
        && iscell(tay.names_of_var) && isempty(tay.names_of_var) ...
        && isnumeric(tay.max_order) && isscalar(tay.max_order) && tay.max_order == 6 ...
        && ischar(tay.opt_method) && strcmp(tay.opt_method,'int') ...
        && isnumeric(tay.eps) && isscalar(tay.eps) && tay.eps == 1e-3 ...
        && isnumeric(tay.tolerance) && isscalar(tay.tolerance) && tay.tolerance == 1e-8;

% ------------------------------ END OF CODE ------------------------------
