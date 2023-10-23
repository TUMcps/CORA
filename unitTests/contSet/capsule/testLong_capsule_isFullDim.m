function res = testLong_capsule_isFullDim
% testLong_capsule_isFullDim - unit test function of isFullDim
%
% Syntax:
%    res = testLong_capsule_isFullDim
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
% Written:       11-March-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;
    
% random number:    0 <= x <= 0.25 ... generator and radius all-zero
%                   0.25 < x <= 0.5 ... generator is all-zero
%                   0.5 < x <= 0.75 ... radius is zero
%                   0.75 < x <= 1 ... generator and radius non-zero
    
nrTests = 1000;
for i=1:nrTests

    % random dimension
    n = randi([2,30]);
    
    % center irrelevant
    c = zeros(n,1);

    randcase = rand(1);
    if randcase <= 0.25
        g = zeros(n,1); r = 0;
        fullDim = false;
    elseif randcase <= 0.5
        g = zeros(n,1); r = rand(1);
        fullDim = true;
    elseif randcase <= 0.75
        g = randn(n,1); r = 0;
        fullDim = false;
    else
        g = randn(n,1); r = rand(1);
        fullDim = true;
    end

    % instantiate capsule
    C = capsule(c,g,r);

    % compute whether full-dimensional
    CfullDim = isFullDim(C);

    if CfullDim ~= fullDim
        res = false;
        break;
    end

end

% ------------------------------ END OF CODE ------------------------------
