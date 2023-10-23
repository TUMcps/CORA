function res = testLong_interval_isFullDim
% testLong_interval_isFullDim - unit test function of isFullDim
%
% Syntax:
%    res = testLong_interval_isFullDim
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

% Authors:       Mark Wetzlinger
% Written:       12-March-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Random cases
res = true;
nrOfTests = 1000;

for i=1:nrOfTests

    % random dimension
    n = randi(50);

    % init random full-dimensional interval
    lb = -rand(n,1);
    ub = rand(n,1);
    I = interval(lb,ub);

    % check with correct solution
    if ~isFullDim(I)
        res = false; break;
    end

    % init random lower-dimensional interval
    % ... by setting random dimension to 0
    randDim = randi(n);
    lb(randDim) = 0;
    ub(randDim) = 0;
    I = interval(lb,ub);

    % check with correct solution
    if isFullDim(I)
        res = false; break;
    end

end

% ------------------------------ END OF CODE ------------------------------
