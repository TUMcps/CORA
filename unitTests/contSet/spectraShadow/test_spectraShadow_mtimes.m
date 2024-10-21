function res = test_spectraShadow_mtimes
% test_spectraShadow_mtimes - unit test function of and
%
% Syntax:
%    res = test_spectraShadow_mtimes
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
 
% assume true
res = true;

A0 = eye(4);
A1 = blkdiag([1 0;0 -1],zeros(2));
A2 = blkdiag(zeros(2),[1 0;0 -1]);

P = rand(2,2);

SpS = spectraShadow([A0,A1,A2],[0;0],P);

M = rand(2,2);

MSpS = M*SpS;

assert(full(all(all(MSpS.G == M*P))));

% ------------------------------ END OF CODE ------------------------------
