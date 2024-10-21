function res = test_spectraShadow_contains
% test_spectraShadow_contains - unit test function of and
%
% Syntax:
%    res = test_spectraShadow_contains
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

% Authors:       Maximilian Perschl, Adrian Kulmburg
% Written:       11-August-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Explicit example
A0 = eye(3);
Ai{1} = [0 1 0;1 0 0;0 0 0];
Ai{2} = [0 0 1;0 0 0;1 0 0];

SpS = spectraShadow([A0,Ai{1},Ai{2}]); %- unit ball

for i = 1:10
    curr_point = -1 + 2*rand(2,1);
    contains_bool = contains(SpS,curr_point,'exact',1e-5);

    assert(contains_bool == (norm(curr_point) <= 1));
end


res = true;

% ------------------------------ END OF CODE ------------------------------
