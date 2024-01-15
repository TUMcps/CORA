function hyp = Inf(n)
% Inf - instantiates a fullspace constrained hyperplane
%
% Syntax:
%    hyp = Inf(n)
%
% Inputs:
%    n - dimension
%
% Outputs:
%    hyp - conHyperplane object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       09-January-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if n <= 0
    throw(CORAerror('CORA:wrongValue','first','positive'));
end

% equality constraint 0*x = 0 fulfilled for all x
hyp = conHyperplane(zeros(1,n),0);

% ------------------------------ END OF CODE ------------------------------
