function res = test_spectraShadow_plus
% test_spectraShadow_plus - unit test function of plus
%
% Syntax:
%    res = test_spectraShadow_plus
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

% Authors:       Maximilian Perschl
% Written:       11-August-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

A0 = eye(4);
Ai{1} = blkdiag([1 0;0 -1],zeros(2));
Ai{2} = blkdiag(zeros(2),[1 0;0 -1]);

SpS1 = spectraShadow([A0,Ai{1},Ai{2}]); %- unit-box

SpS2 = spectraShadow([A0,Ai{1},Ai{2}],[5;5],3*eye(2)); % enlarged box centered at [5;5]


SpS_sum = SpS1 + SpS2; % -box with radius 4 around [5;5]

for i = 1:10
    curr_point = 10*rand(2,1);
    contains_bool = contains(SpS_sum,curr_point,'exact',1e-5);
    should_be_in = (max(curr_point) <= 9) && (min(curr_point) >= 1);

    assert(contains_bool == should_be_in);
end

res = true;

% ------------------------------ END OF CODE ------------------------------
