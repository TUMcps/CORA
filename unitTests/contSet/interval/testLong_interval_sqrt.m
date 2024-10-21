function res = testLong_interval_sqrt
% testLong_interval_sqrt - unit test function of sqrt
%
% Syntax:
%    res = testLong_interval_sqrt
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
% Written:       08-August-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% random tests
nrOfTests = 1000;

for i=1:nrOfTests

    % if inf < 0, result is NaN -> should throw an error
    I = interval(-randi(1000,1,1)*rand(1),randi(1000,1,1)*rand(1));
    assertThrowsAs(@sqrt,'CORA:outOfDomain',I);

    lb = randi(1000,1,1)*rand(1);
    I = interval(lb,lb*2);
    I_sqrt = sqrt(I);

    assertLoop((I_sqrt.inf > 0),i)
end

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
