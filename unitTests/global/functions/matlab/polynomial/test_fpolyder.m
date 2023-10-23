function res = test_fpolyder
% test_fpolyder - unit test function for fpolyder
%
% Syntax:
%    res = test_fpolyder()
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
% See also: fpolyder

% Authors:       Tobias Ladner
% Written:       30-May-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% run random tests
nrTests = 1000;
resvec = false(1,nrTests); % assume false

for i=1:nrTests
    % random dimension
    order = randi(10);
    p = rand(1,order);

    dp1 = fpolyder(p);
    dp2 = polyder(p);

    resvec(i) = all(withinTol(dp1,dp2));
end

% gather results
res = all(resvec);


% ------------------------------ END OF CODE ------------------------------
