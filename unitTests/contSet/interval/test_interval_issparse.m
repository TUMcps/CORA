function res = test_interval_issparse
% test_interval_issparse - unit test function of issparse
%
% Syntax:
%    res = test_interval_issparse
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Tobias Ladner
% Written:       06-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

resvec = true(0);

% test both sparse
I = interval(speye(2), 2*speye(2));
resvec(end+1) = issparse(I);

% test either sparse
I = interval(speye(2), 2*eye(2));
resvec(end+1) = issparse(I);

I = interval(eye(2), 2*speye(2));
resvec(end+1) = issparse(I);

% test not sparse
I = interval([1;2], [3;4]);
resvec(end+1) = ~issparse(I);

res = all(resvec);


% ------------------------------ END OF CODE ------------------------------
