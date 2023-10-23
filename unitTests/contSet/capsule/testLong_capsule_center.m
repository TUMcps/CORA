function res = testLong_capsule_center
% testLong_capsule_center - unit test function of center
%
% Syntax:
%    res = testLong_capsule_center
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
% Written:       28-August-2019
% Last update:   12-March-2021 (add empty case)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------


% random capsules
res = true;
for n=1:50
    
    % init capsule
    c_true = randn(n,1);
    C = capsule(c_true, randn(n,1), rand(1));

    % read center
    c = center(C);
    
    % check result
    if ~compareMatrices(c,c_true)
        res = false; break;
    end
end

% ------------------------------ END OF CODE ------------------------------
