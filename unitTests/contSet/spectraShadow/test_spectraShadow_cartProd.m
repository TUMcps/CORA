function res = test_spectraShadow_cartProd
% test_spectraShadow_cartProd - unit test function of cartProd
%
% Syntax:
%    res = test_spectraShadow_cartProd
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

A0 = eye(3);
Ai{1} = [0 1 0;1 0 0;0 0 0];
Ai{2} = [0 0 1;0 0 0;1 0 0];

SpS1 = spectraShadow([A0,Ai{1},Ai{2}]); %- unit ball

B0 = eye(4);
Bi{1} = blkdiag([1 0;0 -1],zeros(2));
Bi{2} = blkdiag(zeros(2),[1 0;0 -1]);

SpS2 = spectraShadow([B0,Bi{1},Bi{2}]); %- unit-box

SpS_cartProd = cartProd(SpS1,SpS2); 

for i = 1:10
    curr_point_12 = -1 + 2*rand(2,1);
    curr_point_34 = -1 + 2*rand(2,1);
    contains_bool = contains(SpS_cartProd,[curr_point_12;curr_point_34],'exact',1e-5);
    
    should_be_in = (norm(curr_point_12) <= 1) && (max(abs(curr_point_34)) <= 1);

    assert(should_be_in == contains_bool);
end


res = true;

% ------------------------------ END OF CODE ------------------------------
