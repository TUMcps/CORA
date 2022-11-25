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
%    res - boolean 
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Mark Wetzlinger
% Written:      29-August-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

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
res_withResult = Int_reshape == Int_true;

% wrong reshape dimensions should throw error
reshapeSize = [9, 2];
try
    reshape(Int, reshapeSize);
    res_error = false;
catch
    res_error = true;
end

% reshaping to larger array dimension should throw error
reshapeSize = [4, 3, 2];
try
    reshape(Int, reshapeSize);
    res_diffDim = false;
catch
    res_diffDim = true;
end


% add results
res = res_withResult && res_error && res_diffDim;

if res
    disp('test_reshape successful');
else
    disp('test_reshape failed');
end

%------------- END OF CODE --------------