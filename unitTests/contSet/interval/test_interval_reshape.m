function res = test_interval_reshape
% test_interval_reshape - unit test function of reshape
%
% Syntax:
%    res = test_interval_reshape
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
% Written:       29-August-2019
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% create interval
initSize = [6, 4];
lower = -rand(initSize);
upper = rand(initSize);
Int = interval(lower, upper);

% compute reshape
reshapeSize = [8, 3];
Int_reshape = reshape(Int, reshapeSize);

% true solution
lower_true = reshape(lower, reshapeSize);
upper_true = reshape(upper, reshapeSize);
Int_true = interval(lower_true, upper_true);

% compare results
res = Int_reshape == Int_true;

% wrong reshape dimensions should throw error
if CHECKS_ENABLED

reshapeSize = [9, 2];
try
    reshape(Int, reshapeSize);
    res = false;
end

% reshaping to larger array dimension should throw error
reshapeSize = [4, 3, 2];
try
    reshape(Int, reshapeSize);
    res = false;
end

end


% add results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
