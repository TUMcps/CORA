function hyp = empty(n)
% empty - instantiates an empty constrained hyperplane
%
% Syntax:
%    hyp = empty(n)
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
% Last update:   15-January-2024 (TL, parse input)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input
if nargin == 0
    n = 0;
end
inputArgsCheck({{n,'att','numeric',{'scalar','nonnegative'}}});

% equality constraint 0*x = -1 cannot be fulfilled
hyp = conHyperplane(zeros(1,n),-1);

% ------------------------------ END OF CODE ------------------------------
