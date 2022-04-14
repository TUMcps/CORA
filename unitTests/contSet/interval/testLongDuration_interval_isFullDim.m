function res = testLongDuration_interval_isFullDim
% testLongDuration_interval_isFullDim - unit test function of isFullDim
%
% Syntax:  
%    res = testLongDuration_interval_isFullDim
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
% Written:      12-March-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

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


if res
    disp('test_isFullDim successful');
else
    disp('test_isFullDim failed');
end

%------------- END OF CODE --------------
