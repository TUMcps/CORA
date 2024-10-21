function res = testLong_spectraShadow_dim
% testLong_spectraShadow_dim - unit test function of dim
%
% Syntax:
%    res = testLong_spectraShadow_dim
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

% Authors:       Adrian Kulmburg
% Written:       14-August-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

runs = 5;

for i=1:runs
    dimension = randi(8);
    SpS = spectraShadow.generateRandom('Dimension',dimension);
    assert(dim(SpS) == dimension)
end

res = true;

% ------------------------------ END OF CODE ------------------------------
